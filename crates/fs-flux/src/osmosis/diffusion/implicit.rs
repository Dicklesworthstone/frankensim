//! Backward Euler for one solute on the original solvent-fixed molecular face.
//!
//! One conserved transfer q solves q = dt*A*N(nL-q, nR+q). The convex binary
//! free energy and positive pair friction put a physical root between zero and
//! the concentration-equalizing transfer. Material diffusivity may be nonlinear;
//! no uniqueness theorem for arbitrary source curves is asserted. Every trial
//! calls the original constitutive owner at its actual finite inventories.

use super::{
    AqueousSolution, OsmoticNetwork, OsmoticNetworkError, SolutionDiffusionLink,
    SolutionDiffusionSolve, finite, nonnegative, poll, positive, source,
};

/// Nonlinear amount accuracy and work, independent of temporal accuracy and
/// the caller's physical withdrawal/free-energy/conservation limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitBinaryDiffusionControls {
    /// Absolute allowance on EACH original backward-Euler amount equation [mol].
    pub absolute_residual_tolerance_mol: f64,
    /// Relative allowance against the concentration-equalizing transfer, in
    /// [0, 1e-3]. A large common inventory cannot conceal a small exchange.
    pub relative_residual_tolerance: f64,
    /// Maximum candidate-root evaluations per pair, after its initial sample.
    pub maximum_iterations: usize,
    /// Maximum molecular flux evaluations for the COMPLETE requested window,
    /// including initial samples, root trials, rejected steps and both half
    /// sweeps. Native source qualification queries remain additional work.
    pub maximum_evaluations: usize,
}

impl Default for ImplicitBinaryDiffusionControls {
    fn default() -> Self {
        Self {
            absolute_residual_tolerance_mol: 1e-20,
            relative_residual_tolerance: 1e-10,
            maximum_iterations: 96,
            maximum_evaluations: 1_000_000,
        }
    }
}

impl ImplicitBinaryDiffusionControls {
    fn admit(self) -> Result<(), OsmoticNetworkError> {
        if !positive(self.absolute_residual_tolerance_mol)
            || !nonnegative(self.relative_residual_tolerance)
            || self.relative_residual_tolerance > 1e-3
            || self.maximum_iterations == 0
            || self.maximum_iterations > 256
            || self.maximum_evaluations == 0
        {
            return Err(OsmoticNetworkError::Invalid(
                "binary implicit diffusion needs positive absolute amount accuracy, relative accuracy in [0,1e-3], 1..256 root trials and a positive whole-window evaluation budget",
            ));
        }
        Ok(())
    }
}

/// Reusable original thermodynamic inventories, never accepted physical state.
#[derive(Clone, Debug)]
pub(crate) struct ImplicitBinaryDiffusion {
    controls: ImplicitBinaryDiffusionControls,
    trial: [AqueousSolution; 2],
    evaluations: usize,
}

impl OsmoticNetwork {
    /// Select binary backward-Euler molecular stages for an existing fixed-
    /// volume or coupled water/species network. Exactly ONE original solute
    /// is required; the explicit constructors keep their existing behavior.
    ///
    /// Prepares two reusable original solution buffers within the total storage
    /// budget. The same full/two-half timeline controls first-order temporal and
    /// pair-splitting error; implicit stability does not certify time accuracy.
    /// Every accepted transfer still obeys `maximum_species_fraction`: callers
    /// may deliberately relax that physical bound independently of stability.
    /// No nonlinear or source failure publishes physical state or time.
    pub fn with_implicit_binary_diffusion(
        mut self,
        controls: ImplicitBinaryDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        controls.admit()?;
        if self.species.len() != 1
            || self.diffusion_controls.is_none()
            || self.implicit_multicomponent.is_some()
        {
            return Err(OsmoticNetworkError::Invalid(
                "implicit binary diffusion requires one solute and prepared molecular links",
            ));
        }
        if let Some(prepared) = &mut self.implicit_diffusion {
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
        let extra = self.cells[0]
            .solution
            .heap_bytes()
            .checked_mul(2)
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<ImplicitBinaryDiffusion>()))
            .ok_or(OsmoticNetworkError::Invalid(
                "implicit binary solution storage overflow",
            ))?;
        let required =
            self.workspace_bytes
                .checked_add(extra)
                .ok_or(OsmoticNetworkError::Invalid(
                    "implicit binary workspace byte count overflow",
                ))?;
        if required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        self.implicit_diffusion = Some(ImplicitBinaryDiffusion {
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

    /// Deliberately selected molecular integration, distinct from temporal and
    /// physical budgets. None denotes the original explicit molecular stages.
    pub fn implicit_binary_diffusion_controls(&self) -> Option<ImplicitBinaryDiffusionControls> {
        self.implicit_diffusion.as_ref().map(|state| state.controls)
    }
}

impl ImplicitBinaryDiffusion {
    pub(crate) fn begin_window(&mut self) {
        self.evaluations = 0;
    }

    fn sample(
        &mut self,
        link: &mut SolutionDiffusionLink,
        controls: SolutionDiffusionSolve,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, OsmoticNetworkError> {
        poll(cancelled)?;
        if self.evaluations >= self.controls.maximum_evaluations {
            return Err(OsmoticNetworkError::WorkBudget);
        }
        self.evaluations += 1;
        let mut flux = [0.0];
        link.flux_into(
            &self.trial[0],
            &self.trial[1],
            controls,
            &mut flux,
            evaluations,
            cancelled,
        )?;
        Ok(flux[0])
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn transfer(
        &mut self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        link: &mut SolutionDiffusionLink,
        dt: f64,
        controls: SolutionDiffusionSolve,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, OsmoticNetworkError> {
        let amounts = [left.solute_moles()[0], right.solute_moles()[0]];
        let volumes = [left.solvent_volume_m3(), right.solvent_volume_m3()];
        for i in 0..2 {
            self.trial[i]
                .replace_inventories(volumes[i], &amounts[i..i + 1])
                .map_err(source)?;
        }
        let initial_flux = self.sample(link, controls, evaluations, cancelled)?;
        let delta = amounts[0] / volumes[0] - amounts[1] / volumes[1];
        if initial_flux == 0.0 {
            return if delta == 0.0 {
                Ok(0.0)
            } else {
                Err(OsmoticNetworkError::StepUnderflow)
            };
        }
        if !positive(dt)
            || !positive(link.area_m2)
            || !delta.is_finite()
            || delta == 0.0
            || delta.is_sign_positive() != initial_flux.is_sign_positive()
        {
            return Err(OsmoticNetworkError::Invalid(
                "binary molecular flux must follow its original convex concentration force",
            ));
        }
        // V_L V_R / (V_L+V_R), without overflowing the product or sum.
        let (small, large) = (volumes[0].min(volumes[1]), volumes[0].max(volumes[1]));
        let equilibrium = finite((small / (1.0 + small / large)) * delta)?;
        let scale = equilibrium.abs();
        if !positive(scale) {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        let direction = equilibrium.signum();
        let tolerance = finite(
            self.controls.absolute_residual_tolerance_mol
                + self.controls.relative_residual_tolerance * scale,
        )?;
        let initial = integrated_flux(dt, link.area_m2, initial_flux).abs();
        if initial == 0.0 {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        let (mut lower, mut upper) = (0.0, scale);
        // Exact for linear binary diffusion; otherwise only a bounded first
        // trial. Its saturated form never overflows scale + initial.
        let estimate = if initial.is_finite() {
            if initial > scale {
                scale / (1.0 + scale / initial)
            } else {
                initial / (1.0 + initial / scale)
            }
        } else {
            0.5 * scale
        };
        let mut candidate = if estimate > lower && estimate < upper {
            estimate
        } else {
            0.5 * scale
        };
        let mut last_residual = initial;
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            if candidate <= lower || candidate >= upper {
                break;
            }
            let transfer = direction * candidate;
            let next = [
                finite(amounts[0] - transfer)?,
                finite(amounts[1] + transfer)?,
            ];
            if next.iter().any(|&value| !nonnegative(value)) {
                return Err(OsmoticNetworkError::Invalid(
                    "binary equilibrium transfer leaves the nonnegative inventory bracket",
                ));
            }
            for i in 0..2 {
                self.trial[i]
                    .replace_inventories(volumes[i], &next[i..i + 1])
                    .map_err(source)?;
            }
            let flux = self.sample(link, controls, evaluations, cancelled)?;
            let predicted = integrated_flux(dt, link.area_m2, flux);
            // Check BOTH actual rounded changes, not merely the nominal root
            // coordinate. Infinite exposure or a zero rounded flux at equilibrium
            // cannot create an Inf<=Inf or an underflowed rate-form false root.
            let residuals = [
                (amounts[0] - next[0]) - predicted,
                (next[1] - amounts[1]) - predicted,
            ];
            last_residual = residuals[0].abs().max(residuals[1].abs());
            if predicted.is_finite()
                && residuals.iter().all(|v| v.is_finite())
                && last_residual <= tolerance
            {
                if next == amounts {
                    return Err(OsmoticNetworkError::StepUnderflow);
                }
                poll(cancelled)?;
                return Ok(transfer);
            }
            let signed_residual = (transfer - predicted) * direction;
            if signed_residual > 0.0 {
                upper = candidate;
            } else {
                lower = candidate;
            }
            candidate = lower + 0.5 * (upper - lower);
        }
        Err(OsmoticNetworkError::ImplicitDiffusionConvergence {
            residual_mol: last_residual,
            maximum_mol: tolerance,
        })
    }
}

// Multiply three finite factors by pairing the smallest and largest magnitudes
// first. A finite representable integral is not lost to dt*A overflow or an
// A*N underflow. A genuinely overflowing prediction may bracket a root but is
// NEVER accepted as a residual. Exact zero flux has exact zero integral.
pub(super) fn integrated_flux(dt: f64, area: f64, flux: f64) -> f64 {
    if flux == 0.0 {
        return 0.0;
    }
    let mut factors = [dt, area, flux.abs()];
    factors.sort_by(f64::total_cmp);
    ((factors[0] * factors[2]) * factors[1]).copysign(flux)
}
