//! Conservative multisolute backward Euler through the native Newton owner.
//!
//! Independent transfers retain all original Maxwell--Stefan cross terms.
//! Bounded coordinates keep trial inventories nonnegative; only the physical
//! owner and final molar residual gates can admit an endpoint. This adapter
//! supplies no competing Newton, Krylov or thermodynamic implementation.

use super::implicit::integrated_flux;
use super::{
    AqueousSolution, MAX_OSMOTIC_SPECIES, OsmoticNetwork, OsmoticNetworkError,
    SolutionDiffusionLink, SolutionDiffusionSolve, finite, nonnegative, poll, positive, source,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonKrylovConfig, NewtonKrylovState, NewtonStallDiagnosis, SolverRunProgress};
use std::cell::RefCell;

/// Physical molar accuracy and bounded native nonlinear work for one complete
/// requested window. Temporal error and physical withdrawal remain separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitMulticomponentDiffusionControls {
    /// Absolute allowance on EACH actual endpoint amount equation [mol].
    pub absolute_residual_tolerance_mol: f64,
    /// Additional allowance times the actual transfer magnitude, in [0,1e-3].
    /// Common background inventory never enlarges this allowance.
    pub relative_residual_tolerance: f64,
    /// Native Newton attempts per physical pair, in 1..=64.
    pub maximum_iterations: usize,
    /// All original molecular evaluations in the complete requested window:
    /// initial samples, Jacobian probes, line searches, final replay and retries.
    pub maximum_evaluations: usize,
    /// Numerical directional-difference step in bounded logit coordinates.
    pub difference_step: f64,
    /// Native policy on residuals divided by the FIXED absolute molar allowance.
    /// Relative convergence is disabled; the absolute target is at most 0.5.
    /// Physical abs+relative acceptance is checked at every native checkpoint.
    pub newton: NewtonKrylovConfig,
}

impl Default for ImplicitMulticomponentDiffusionControls {
    fn default() -> Self {
        Self {
            absolute_residual_tolerance_mol: 1e-20,
            relative_residual_tolerance: 1e-10,
            maximum_iterations: 64,
            maximum_evaluations: 1_000_000,
            difference_step: 1e-5,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 0.25,
                relative_tolerance: 0.0,
                linear_restart: MAX_OSMOTIC_SPECIES,
                max_linear_cycles: 4,
                // Directional finite differences are numerical actions, so an
                // inner solve demanding analytic-tangent precision can stall.
                forcing_minimum: 1e-3,
                ..NewtonKrylovConfig::default()
            },
        }
    }
}

impl ImplicitMulticomponentDiffusionControls {
    /// Conservative payload allowance for the existing allocating Newton/FD/
    /// FGMRES implementation: at most 64 unknowns, restart 64, 8 linear cycles
    /// and 64 outer attempts, including staged clones and telemetry. This is
    /// an admission estimate, NOT an allocated arena or process-memory bound.
    /// Allocator metadata and transient original source queries are additional.
    pub const fn solver_workspace_bytes() -> usize {
        512 * 1024
    }

    fn admit(self) -> Result<(), OsmoticNetworkError> {
        self.newton.validate().map_err(|_| {
            OsmoticNetworkError::Invalid("invalid native multicomponent Newton controls")
        })?;
        if !positive(self.absolute_residual_tolerance_mol)
            || !nonnegative(self.relative_residual_tolerance)
            || self.relative_residual_tolerance > 1e-3
            || !(1..=64).contains(&self.maximum_iterations)
            || self.maximum_evaluations == 0
            || !positive(self.difference_step)
            || self.difference_step > 0.1
            || self.newton.absolute_tolerance > 0.5
            || self.newton.relative_tolerance != 0.0
            || self.newton.linear_restart > MAX_OSMOTIC_SPECIES
            || self.newton.max_linear_cycles > 8
        {
            return Err(OsmoticNetworkError::Invalid(
                "invalid multicomponent molar accuracy, difference scale or bounded nonlinear work",
            ));
        }
        Ok(())
    }

    fn tolerance(self, transfer: f64) -> Result<f64, OsmoticNetworkError> {
        finite(
            self.absolute_residual_tolerance_mol
                + self.relative_residual_tolerance * transfer.abs(),
        )
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ImplicitMulticomponentDiffusion {
    controls: ImplicitMulticomponentDiffusionControls,
    trial: [AqueousSolution; 2],
    evaluations: usize,
}

impl OsmoticNetwork {
    /// Select fully coupled backward-Euler molecular pair equations with two
    /// or more original solutes. Applies to fixed-volume and coupled water/
    /// species networks on their existing adaptive timeline. The specialized
    /// one-solute builder and explicit constructors retain their behavior.
    ///
    /// Original friction, chemical forces, pinned sources and every physical
    /// acceptance gate remain authoritative. Native solver work allocates
    /// transient vectors inside the declared conservative payload allowance.
    pub fn with_implicit_multicomponent_diffusion(
        mut self,
        controls: ImplicitMulticomponentDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        controls.admit()?;
        if self.species.len() < 2
            || self.diffusion_controls.is_none()
            || self.implicit_diffusion.is_some()
            || self.zero_current_diffusion().is_some()
            || self.implicit_salt.is_some()
            || self.electroneutral_diffusion.is_some()
        {
            return Err(OsmoticNetworkError::Invalid(
                "multicomponent implicit diffusion needs at least two original neutral solutes and prepared molecular links",
            ));
        }
        if let Some(prepared) = &mut self.implicit_multicomponent {
            if self.workspace_bytes > maximum_workspace_bytes {
                return Err(OsmoticNetworkError::Workspace {
                    required_bytes: self.workspace_bytes,
                    budget_bytes: maximum_workspace_bytes,
                });
            }
            prepared.controls = controls;
            prepared.begin_window();
            return Ok(self);
        }
        let required = self.cells[0]
            .solution
            .heap_bytes()
            .checked_mul(2)
            .and_then(|n| n.checked_add(core::mem::size_of::<ImplicitMulticomponentDiffusion>()))
            .and_then(|n| {
                n.checked_add(ImplicitMulticomponentDiffusionControls::solver_workspace_bytes())
            })
            .and_then(|n| n.checked_add(self.workspace_bytes))
            .ok_or(OsmoticNetworkError::Invalid(
                "multicomponent implicit workspace overflow",
            ))?;
        if required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        self.implicit_multicomponent = Some(ImplicitMulticomponentDiffusion {
            controls,
            trial: [
                self.cells[0].solution.clone(),
                self.cells[0].solution.clone(),
            ],
            evaluations: 0,
        });
        self.workspace_bytes = required;
        Ok(self)
    }

    pub fn implicit_multicomponent_diffusion_controls(
        &self,
    ) -> Option<ImplicitMulticomponentDiffusionControls> {
        self.implicit_multicomponent
            .as_ref()
            .map(|state| state.controls)
    }
}

struct PairLayout {
    count: usize,
    volumes: [f64; 2],
    amounts: [[f64; MAX_OSMOTIC_SPECIES]; 2],
    bounds: [[f64; 2]; MAX_OSMOTIC_SPECIES],
    active: Vec<usize>,
}

struct Candidate {
    transfers: [f64; MAX_OSMOTIC_SPECIES],
    residuals: [[f64; MAX_OSMOTIC_SPECIES]; 2],
    scaled: Vec<f64>,
    changed: bool,
}

impl ImplicitMulticomponentDiffusion {
    pub(crate) fn begin_window(&mut self) {
        self.evaluations = 0;
    }

    fn sample(
        &mut self,
        link: &mut SolutionDiffusionLink,
        controls: SolutionDiffusionSolve,
        evaluations: &mut usize,
        flux: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        poll(cancelled)?;
        if self.evaluations >= self.controls.maximum_evaluations {
            return Err(OsmoticNetworkError::WorkBudget);
        }
        self.evaluations += 1;
        link.flux_into(
            &self.trial[0],
            &self.trial[1],
            controls,
            flux,
            evaluations,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &mut self,
        layout: &PairLayout,
        x: &[f64],
        link: &mut SolutionDiffusionLink,
        dt: f64,
        controls: SolutionDiffusionSolve,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Candidate, OsmoticNetworkError> {
        poll(cancelled)?;
        let n = layout.count;
        let mut transfers = [0.0; MAX_OSMOTIC_SPECIES];
        for (&coordinate, &i) in x.iter().zip(&layout.active) {
            // Evaluate both logistic weights independently: tiny invading
            // fractions are not erased by subtracting a rounded value from 1.
            let small = (-coordinate.abs()).exp();
            let denominator = 1.0 + small;
            let (left, right) = if coordinate >= 0.0 {
                (small / denominator, 1.0 / denominator)
            } else {
                (1.0 / denominator, small / denominator)
            };
            let [lower, upper] = layout.bounds[i];
            transfers[i] = finite(lower * left + upper * right)?;
        }
        let mut next = layout.amounts;
        let mut changed = false;
        for i in 0..n {
            next[0][i] = finite(layout.amounts[0][i] - transfers[i])?;
            next[1][i] = finite(layout.amounts[1][i] + transfers[i])?;
            if !nonnegative(next[0][i]) || !nonnegative(next[1][i]) {
                return Err(OsmoticNetworkError::Invalid(
                    "bounded molecular trial left nonnegative inventories",
                ));
            }
            changed |= next[0][i] != layout.amounts[0][i] || next[1][i] != layout.amounts[1][i];
        }
        for side in 0..2 {
            self.trial[side]
                .replace_inventories(layout.volumes[side], &next[side][..n])
                .map_err(source)?;
        }
        let mut flux = [0.0; MAX_OSMOTIC_SPECIES];
        self.sample(link, controls, evaluations, &mut flux[..n], cancelled)?;
        let mut residuals = [[0.0; MAX_OSMOTIC_SPECIES]; 2];
        let mut nominal = [0.0; MAX_OSMOTIC_SPECIES];
        for i in 0..n {
            let integrated = finite(integrated_flux(dt, link.area_m2, flux[i]))?;
            if flux[i] != 0.0 && integrated == 0.0 {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            residuals[0][i] = finite((layout.amounts[0][i] - next[0][i]) - integrated)?;
            residuals[1][i] = finite((next[1][i] - layout.amounts[1][i]) - integrated)?;
            // Fixed scaling is essential: dividing by |q| here can turn a
            // weak-step equation into an almost constant sign away from q=0.
            // The physical relative allowance is applied only to checkpoints.
            nominal[i] = finite(
                (transfers[i] - integrated) / self.controls.absolute_residual_tolerance_mol,
            )?;
        }
        Ok(Candidate {
            transfers,
            residuals,
            scaled: layout.active.iter().map(|&i| nominal[i]).collect(),
            changed,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn transfers_into(
        &mut self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        link: &mut SolutionDiffusionLink,
        dt: f64,
        controls: SolutionDiffusionSolve,
        evaluations: &mut usize,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        let numerical = self.controls;
        let n = left.solute_moles().len();
        if output.len() != n || !positive(dt) || !positive(link.area_m2) {
            return Err(OsmoticNetworkError::Invalid(
                "complete multicomponent output and finite positive duration/area required",
            ));
        }
        let mut layout = PairLayout {
            count: n,
            volumes: [left.solvent_volume_m3(), right.solvent_volume_m3()],
            amounts: [[0.0; MAX_OSMOTIC_SPECIES]; 2],
            bounds: [[0.0; 2]; MAX_OSMOTIC_SPECIES],
            active: Vec::with_capacity(n),
        };
        for (side, original) in [left, right].into_iter().enumerate() {
            layout.amounts[side][..n].copy_from_slice(original.solute_moles());
            self.trial[side]
                .replace_inventories(layout.volumes[side], original.solute_moles())
                .map_err(source)?;
        }
        let mut initial_flux = [0.0; MAX_OSMOTIC_SPECIES];
        self.sample(
            link,
            controls,
            evaluations,
            &mut initial_flux[..n],
            cancelled,
        )?;
        if initial_flux[..n].iter().all(|&v| v == 0.0) {
            if (0..n).any(|i| {
                layout.amounts[0][i] / layout.volumes[0] != layout.amounts[1][i] / layout.volumes[1]
            }) {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            poll(cancelled)?;
            output.fill(0.0);
            return Ok(());
        }
        let mut guess = Vec::with_capacity(n);
        let caps = link
            .material()
            .thermodynamics()
            .maximum_concentrations_mol_m3();
        let small_volume = layout.volumes[0].min(layout.volumes[1]);
        let large_volume = layout.volumes[0].max(layout.volumes[1]);
        let effective_volume = small_volume / (1.0 + small_volume / large_volume);
        for (i, &cap) in caps.iter().enumerate() {
            let (a, b) = (layout.amounts[0][i], layout.amounts[1][i]);
            // A product overflowing to infinity simply cannot restrict the
            // finite conserved inventory interval. Include the already-admitted
            // original q=0 against a one-ULP capacity-product discrepancy.
            let lower = (-b).max(a - cap * layout.volumes[0]).min(0.0);
            let upper = a.min(cap * layout.volumes[1] - b).max(0.0);
            layout.bounds[i] = [lower, upper];
            if lower == upper {
                continue;
            }
            let equilibrium = effective_volume * (a / layout.volumes[0] - b / layout.volumes[1]);
            let start = if equilibrium > lower && equilibrium < upper {
                equilibrium
            } else {
                0.5 * lower + 0.5 * upper
            };
            let coordinate = finite((start - lower).ln() - (upper - start).ln())?;
            layout.active.push(i);
            guess.push(coordinate);
        }
        if layout.active.is_empty() {
            return Err(OsmoticNetworkError::Invalid(
                "nonzero molecular force has no concentration-supported transfer interval",
            ));
        }
        let owner = RefCell::new((self, link, evaluations));
        let cancel = RefCell::new(cancelled);
        let evaluate = |x: &[f64]| {
            let mut owner = owner.borrow_mut();
            let (state, link, evaluations) = &mut *owner;
            state.evaluate(&layout, x, link, dt, controls, evaluations, &mut || {
                (*cancel.borrow_mut())()
            })
        };
        let problem = FallibleFiniteDifference::new(
            layout.active.len(),
            numerical.difference_step,
            |x: &[f64]| evaluate(x).map(|candidate| candidate.scaled),
        )
        .map_err(numerical_error)?;
        let construction = NewtonKrylovState::new(&problem, guess, numerical.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical_error(error));
        }
        let mut solver = construction.map_err(|_| {
            OsmoticNetworkError::Invalid("native multicomponent Newton initialization refused")
        })?;
        let mut terminal = None;
        for attempt in 0..=numerical.maximum_iterations {
            // Fresh original-law replays of native checkpoints, never Jacobian
            // probes, supply BOTH actual rounded molar equations. Physical
            // relative accuracy may stop before the stricter absolute target.
            let candidate = evaluate(&solver.x)?;
            let mut physical_refusal = None;
            for i in 0..n {
                let maximum = numerical.tolerance(candidate.transfers[i])?;
                let residual = candidate.residuals[0][i]
                    .abs()
                    .max(candidate.residuals[1][i].abs());
                if !residual.is_finite() || residual > maximum {
                    physical_refusal = Some(OsmoticNetworkError::ImplicitDiffusionConvergence {
                        residual_mol: residual,
                        maximum_mol: maximum,
                    });
                    break;
                }
            }
            if physical_refusal.is_none() {
                if !candidate.changed {
                    return Err(OsmoticNetworkError::StepUnderflow);
                }
                poll(&mut || (*cancel.borrow_mut())())?;
                output.copy_from_slice(&candidate.transfers[..n]);
                return Ok(());
            }
            if terminal == Some(None) {
                return Err(physical_refusal.expect("failed original molar gate"));
            }
            if terminal.is_some() || attempt == numerical.maximum_iterations {
                return Err(OsmoticNetworkError::ImplicitMulticomponentConvergence {
                    scaled_residual: solver.residual_norm(),
                    diagnosis: terminal
                        .flatten()
                        .or(Some(NewtonStallDiagnosis::BudgetExhausted)),
                });
            }
            let result = solver.run_controlled(&problem, 1, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
            if let Some(error) = problem.take_failure() {
                return Err(numerical_error(error));
            }
            match result.progress {
                SolverRunProgress::Complete => {}
                SolverRunProgress::Paused => return Err(OsmoticNetworkError::Cancelled),
                _ => {
                    return Err(OsmoticNetworkError::Invalid(
                        "native multicomponent nonlinear callback interrupted",
                    ));
                }
            }
            if result.report.converged
                || result.report.diagnosis != Some(NewtonStallDiagnosis::BudgetExhausted)
            {
                terminal = Some(result.report.diagnosis);
            }
        }
        unreachable!("bounded native checkpoint loop returns at its final attempt")
    }
}

fn numerical_error(error: FiniteDifferenceError<OsmoticNetworkError>) -> OsmoticNetworkError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(_) => OsmoticNetworkError::Invalid(
            "native multicomponent numerical Jacobian refused its finite stencil",
        ),
    }
}
