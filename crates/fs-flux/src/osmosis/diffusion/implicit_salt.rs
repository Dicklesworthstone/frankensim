//! Backward Euler on the one charge-neutral extent of a declared binary salt.
//!
//! For opposite valences z, the gcd-reduced positive stoichiometry nu satisfies
//! z.nu=0. A single transfer q changes both original ion inventories by nu*q;
//! no nonlinear probe independently perturbs an ion or invents a counterion.
//! Every probe evaluates the original zero-current Maxwell--Stefan owner at
//! its actual rounded endpoint inventories and re-queries original sources.
//! Arbitrary nonlinear source curves do not inherit a uniqueness assertion.

use super::implicit::{TransferResidual, integrated_flux, solve_bounded_transfer};
use super::{
    AqueousSolution, OsmoticNetwork, OsmoticNetworkError, SolutionDiffusionLink,
    SolutionDiffusionSolve, finite, nonnegative, poll, positive, source,
};

/// Amount accuracy and bounded work for a two-ion charge-neutral exchange.
/// Temporal error and physical withdrawal limits retain their separate owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitSaltDiffusionControls {
    /// Absolute allowance on EACH original ion's endpoint amount equation [mol].
    pub absolute_residual_tolerance_mol: f64,
    /// Relative allowance against that ion's actual proposed transfer, in
    /// [0,1e-3]. Background salt cannot conceal a small exchange residual.
    pub relative_residual_tolerance: f64,
    /// Maximum candidate evaluations per pair, including bracket qualification
    /// after the initial sample. Bounds are 1..256; no concentration floor is used.
    pub maximum_iterations: usize,
    /// Maximum original molecular evaluations for the COMPLETE requested window,
    /// including rejected steps and full/two-half comparisons. Original source
    /// qualification queries remain additional work.
    pub maximum_evaluations: usize,
}

impl Default for ImplicitSaltDiffusionControls {
    fn default() -> Self {
        Self {
            absolute_residual_tolerance_mol: 1e-20,
            relative_residual_tolerance: 1e-10,
            maximum_iterations: 96,
            maximum_evaluations: 1_000_000,
        }
    }
}

impl ImplicitSaltDiffusionControls {
    fn admit(self) -> Result<(), OsmoticNetworkError> {
        if !positive(self.absolute_residual_tolerance_mol)
            || !nonnegative(self.relative_residual_tolerance)
            || self.relative_residual_tolerance > 1e-3
            || !(1..=256).contains(&self.maximum_iterations)
            || self.maximum_evaluations == 0
        {
            return Err(OsmoticNetworkError::Invalid(
                "implicit salt diffusion needs positive absolute ion accuracy, relative accuracy in [0,1e-3], 1..256 root trials and a positive whole-window evaluation budget",
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
pub(crate) struct ImplicitSaltDiffusion {
    controls: ImplicitSaltDiffusionControls,
    stoichiometry: [f64; 2],
    trial: [AqueousSolution; 2],
    evaluations: usize,
}

impl OsmoticNetwork {
    /// Select charge-constrained backward Euler for exactly two declared ions
    /// with opposite nonzero valences on an existing zero-current forest.
    /// Their gcd-reduced valences define the ONE salt extent; no dissociation,
    /// neutral spectator, independent salt or electrode process is inferred.
    ///
    /// Applies to fixed-volume diffusion and fully retaining water/molecular
    /// transport on the original adaptive full/two-half timeline. That timeline
    /// still controls temporal and pair-splitting error. Every accepted ion
    /// transfer keeps the caller's withdrawal, charge, conservation, source and
    /// free-energy gates. Implicit stability does not waive those limits.
    ///
    /// Two reusable original solution buffers are included in the total storage
    /// budget. Constant-data pair work allocates no vectors. A source, nonlinear
    /// or cancellation failure publishes neither physical state nor time.
    pub fn with_implicit_salt_diffusion(
        mut self,
        controls: ImplicitSaltDiffusionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        controls.admit()?;
        if self.species.len() != 2
            || self.diffusion_controls.is_none()
            || self.implicit_diffusion.is_some()
            || self.implicit_multicomponent.is_some()
            || self.electroneutral_diffusion.is_some()
        {
            return Err(OsmoticNetworkError::Invalid(
                "implicit salt diffusion requires exactly two ions on a prepared zero-current forest",
            ));
        }
        let valences = self
            .zero_current_diffusion()
            .ok_or(OsmoticNetworkError::Invalid(
                "implicit salt diffusion requires an explicit zero-current ionic assignment",
            ))?
            .valences();
        if valences.len() != 2
            || valences[0] == 0
            || valences[1] == 0
            || valences[0].is_positive() == valences[1].is_positive()
        {
            return Err(OsmoticNetworkError::Invalid(
                "implicit salt diffusion requires exactly two opposite nonzero valences",
            ));
        }
        // Promote before abs so i32::MIN is an admitted explicit valence too.
        let magnitudes = [
            i64::from(valences[0]).unsigned_abs(),
            i64::from(valences[1]).unsigned_abs(),
        ];
        let (mut a, mut b) = (magnitudes[0], magnitudes[1]);
        while b != 0 {
            (a, b) = (b, a % b);
        }
        let stoichiometry = [(magnitudes[1] / a) as f64, (magnitudes[0] / a) as f64];
        if let Some(prepared) = &mut self.implicit_salt {
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
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<ImplicitSaltDiffusion>()))
            .and_then(|bytes| bytes.checked_add(self.workspace_bytes))
            .ok_or(OsmoticNetworkError::Invalid(
                "implicit salt workspace byte count overflow",
            ))?;
        if required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        self.implicit_salt = Some(ImplicitSaltDiffusion {
            controls,
            stoichiometry,
            trial: [
                self.cells[0].solution.clone(),
                self.cells[0].solution.clone(),
            ],
            evaluations: 0,
        });
        self.workspace_bytes = required;
        Ok(self)
    }

    pub fn implicit_salt_diffusion_controls(&self) -> Option<ImplicitSaltDiffusionControls> {
        self.implicit_salt.as_ref().map(|state| state.controls)
    }
}

struct SaltPair {
    amounts: [[f64; 2]; 2],
    volumes: [f64; 2],
    direction: f64,
    root_ion: usize,
}

impl ImplicitSaltDiffusion {
    pub(crate) fn begin_window(&mut self) {
        self.evaluations = 0;
    }

    fn sample(
        &mut self,
        link: &mut SolutionDiffusionLink,
        controls: SolutionDiffusionSolve,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], OsmoticNetworkError> {
        poll(cancelled)?;
        if self.evaluations >= self.controls.maximum_evaluations {
            return Err(OsmoticNetworkError::WorkBudget);
        }
        self.evaluations += 1;
        let mut flux = [0.0; 2];
        link.flux_into(
            &self.trial[0],
            &self.trial[1],
            controls,
            &mut flux,
            evaluations,
            cancelled,
        )?;
        Ok(flux)
    }

    #[allow(clippy::too_many_arguments)]
    fn residual(
        &mut self,
        pair: &SaltPair,
        coordinate: f64,
        link: &mut SolutionDiffusionLink,
        dt: f64,
        controls: SolutionDiffusionSolve,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransferResidual, OsmoticNetworkError> {
        poll(cancelled)?;
        let transfers = self
            .stoichiometry
            .map(|nu| nu * pair.direction * coordinate);
        let mut next = pair.amounts;
        for i in 0..2 {
            next[0][i] = finite(pair.amounts[0][i] - transfers[i])?;
            next[1][i] = finite(pair.amounts[1][i] + transfers[i])?;
            if !nonnegative(next[0][i]) || !nonnegative(next[1][i]) {
                return Err(OsmoticNetworkError::Invalid(
                    "implicit salt transfer leaves the nonnegative stoichiometric bracket",
                ));
            }
        }
        for side in 0..2 {
            self.trial[side]
                .replace_inventories(pair.volumes[side], &next[side])
                .map_err(source)?;
        }
        // The original owner admits BOTH actual rounded endpoint charges before
        // any ionic face/source query. No projection of its flux is published.
        let flux = self.sample(link, controls, evaluations, cancelled)?;
        let predicted = flux.map(|n| integrated_flux(dt, link.area_m2, n));
        let mut worst = TransferResidual {
            signed_residual: (transfers[pair.root_ion] - predicted[pair.root_ion]) * pair.direction,
            residual_mol: 0.0,
            maximum_mol: self.controls.absolute_residual_tolerance_mol,
        };
        for i in 0..2 {
            let maximum = self.controls.tolerance(transfers[i])?;
            let residuals = [
                (pair.amounts[0][i] - next[0][i]) - predicted[i],
                (next[1][i] - pair.amounts[1][i]) - predicted[i],
            ];
            let residual =
                if predicted[i].is_finite() && residuals.iter().all(|value| value.is_finite()) {
                    residuals[0].abs().max(residuals[1].abs())
                } else {
                    f64::INFINITY
                };
            // Select the actual equation with the largest normalized defect,
            // preserving its own molar tolerance for an honest failure report.
            if residual / maximum >= worst.residual_mol / worst.maximum_mol {
                worst.residual_mol = residual;
                worst.maximum_mol = maximum;
            }
        }
        if worst.converged()
            && (0..2).any(|i| {
                transfers[i] != 0.0
                    && next[0][i] == pair.amounts[0][i]
                    && next[1][i] == pair.amounts[1][i]
            })
        {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        Ok(worst)
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
        cancelled: &mut impl FnMut() -> bool,
        output: &mut [f64],
    ) -> Result<(), OsmoticNetworkError> {
        if output.len() != 2
            || left.solute_moles().len() != 2
            || right.solute_moles().len() != 2
            || link.zero_current_diffusion().is_none()
            || !positive(dt)
            || !positive(link.area_m2)
        {
            return Err(OsmoticNetworkError::Invalid(
                "implicit salt pair needs its original two-ion assignment, positive area and time",
            ));
        }
        let mut pair = SaltPair {
            amounts: [
                [left.solute_moles()[0], left.solute_moles()[1]],
                [right.solute_moles()[0], right.solute_moles()[1]],
            ],
            volumes: [left.solvent_volume_m3(), right.solvent_volume_m3()],
            direction: 1.0,
            root_ion: usize::from(self.stoichiometry[1] < self.stoichiometry[0]),
        };
        for side in 0..2 {
            self.trial[side]
                .replace_inventories(pair.volumes[side], &pair.amounts[side])
                .map_err(source)?;
        }
        let initial_flux = self.sample(link, controls, evaluations, cancelled)?;
        if initial_flux == [0.0; 2] {
            poll(cancelled)?;
            output.fill(0.0);
            return Ok(());
        }
        let initial = integrated_flux(dt, link.area_m2, initial_flux[pair.root_ion]).abs()
            / self.stoichiometry[pair.root_ion];
        if initial == 0.0 {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        pair.direction = initial_flux[pair.root_ion].signum();
        let donor = usize::from(pair.direction < 0.0);
        let limit = (pair.amounts[donor][0] / self.stoichiometry[0])
            .min(pair.amounts[donor][1] / self.stoichiometry[1]);
        if !positive(limit) {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        let (small, large) = (
            pair.volumes[0].min(pair.volumes[1]),
            pair.volumes[0].max(pair.volumes[1]),
        );
        let effective_volume = small / (1.0 + small / large);
        let mut hint: f64 = 0.0;
        for i in 0..2 {
            let delta = finite(
                pair.amounts[0][i] / pair.volumes[0] - pair.amounts[1][i] / pair.volumes[1],
            )?;
            hint = hint.max(finite(
                effective_volume * delta * pair.direction / self.stoichiometry[i],
            )?);
        }
        let (mut lower, mut upper) = (
            0.0,
            if hint > 0.0 && hint < limit {
                hint
            } else {
                0.5 * limit
            },
        );
        let mut remaining = self.controls.maximum_iterations;
        let mut last = TransferResidual {
            signed_residual: -initial,
            residual_mol: initial,
            maximum_mol: self.controls.absolute_residual_tolerance_mol,
        };
        // Qualify the concentration-equalizing bracket with the ORIGINAL law.
        // Tiny admitted endpoint charge residuals can move its equilibrium;
        // bounded interior expansion handles that case without a zero-ion probe
        // or a concentration floor. Arbitrary nonlinear source domains remain
        // authoritative; a refused source query is never retried with stale data.
        let coordinate = loop {
            if remaining == 0 || upper <= lower || upper >= limit {
                return Err(OsmoticNetworkError::ImplicitDiffusionConvergence {
                    residual_mol: last.residual_mol,
                    maximum_mol: last.maximum_mol,
                });
            }
            remaining -= 1;
            last = self.residual(&pair, upper, link, dt, controls, evaluations, cancelled)?;
            if last.converged() {
                break upper;
            }
            if last.signed_residual > 0.0 {
                let estimate = if initial.is_finite() {
                    if initial > upper {
                        upper / (1.0 + upper / initial)
                    } else {
                        initial / (1.0 + initial / upper)
                    }
                } else {
                    lower + 0.5 * (upper - lower)
                };
                break solve_bounded_transfer(
                    [lower, upper],
                    estimate,
                    remaining,
                    last,
                    |candidate| {
                        self.residual(&pair, candidate, link, dt, controls, evaluations, cancelled)
                    },
                )?;
            }
            lower = upper;
            upper += 0.5 * (limit - upper);
        };
        let transfers = self
            .stoichiometry
            .map(|nu| nu * pair.direction * coordinate);
        poll(cancelled)?;
        output.copy_from_slice(&transfers);
        Ok(())
    }
}
