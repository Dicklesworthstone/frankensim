//! Dispatch only: the statistical owners retain sampler, checkpoint and
//! failure semantics. Both paths evaluate the same native child-run callback.

use std::fmt::Display;

use fs_blake3::ContentHash;
use fs_uq::{GaussianCopulaExecution, GaussianCopulaQmcExecution,
    QmcConfig, QmcExecution, QmcReport, SobolExecution, SobolReport, UqExecution, UqStatus};

use super::{Model, Result, compliance, fail, plan};

pub(super) enum Execution {
    SobolSensitivity(SobolExecution),
    MonteCarlo(UqExecution),
    QuasiMonteCarlo(QmcExecution),
    CopulaMonteCarlo(GaussianCopulaExecution),
    CopulaQuasiMonteCarlo(GaussianCopulaQmcExecution),
}

fn layout(model: &Model) -> Option<QmcConfig> {
    model.bound.study().qmc().map(|layout| QmcConfig {
        replicates: layout.replicates,
        samples_per_replicate: layout.samples_per_replicate,
    })
}

impl Execution {
    pub(super) fn new(model: &Model) -> Result<Self> {
        let plan = plan(model);
        if model.bound.study().sobol_sensitivity() {
            return SobolExecution::new(&plan).map(Self::SobolSensitivity)
                .map_err(|error| fail("cli-uncertainty-plan", error));
        }
        match (model.bound.study().latent_correlation(), layout(model)) {
            (Some(matrix), Some(layout)) => GaussianCopulaQmcExecution::new(&plan, matrix, layout)
                .map(Self::CopulaQuasiMonteCarlo),
            (Some(matrix), None) => GaussianCopulaExecution::new(&plan, matrix)
                .map(Self::CopulaMonteCarlo),
            (None, Some(layout)) => QmcExecution::new(&plan, layout).map(Self::QuasiMonteCarlo),
            (None, None) => UqExecution::new(&plan).map(Self::MonteCarlo),
        }.map_err(|error| fail("cli-uncertainty-plan", error))
    }

    pub(super) fn restore(model: &Model, bytes: &[u8]) -> Result<Self> {
        let plan = plan(model);
        if model.bound.study().sobol_sensitivity() {
            return SobolExecution::restore(&plan, model.identity(), bytes).map(Self::SobolSensitivity)
                .map_err(|error| fail("cli-uncertainty-resume", error.to_string()));
        }
        match (model.bound.study().latent_correlation(), layout(model)) {
            (Some(matrix), Some(layout)) => GaussianCopulaQmcExecution::restore(
                &plan, matrix, layout, model.identity(), bytes).map(Self::CopulaQuasiMonteCarlo),
            (Some(matrix), None) => GaussianCopulaExecution::restore(
                &plan, matrix, model.identity(), bytes).map(Self::CopulaMonteCarlo),
            (None, Some(layout)) => QmcExecution::restore(&plan, layout, model.identity(), bytes)
                .map(Self::QuasiMonteCarlo),
            (None, None) => UqExecution::restore(&plan, model.identity(), bytes).map(Self::MonteCarlo),
        }.map_err(|error| fail("cli-uncertainty-resume", error.to_string()))
    }

    pub(super) fn checkpoint(&self, identity: ContentHash) -> Result<Vec<u8>> {
        let bytes = match self {
            Self::SobolSensitivity(execution) => execution.checkpoint(identity),
            Self::MonteCarlo(execution) => execution.checkpoint(identity),
            Self::QuasiMonteCarlo(execution) => execution.checkpoint(identity),
            Self::CopulaMonteCarlo(execution) => execution.checkpoint(identity),
            Self::CopulaQuasiMonteCarlo(execution) => execution.checkpoint(identity),
        };
        bytes.map_err(|error| fail("cli-uncertainty-checkpoint", error.to_string()))
    }

    pub(super) fn observations(&self) -> &[f64] {
        match self {
            Self::SobolSensitivity(execution) => execution.observations(),
            Self::MonteCarlo(execution) => execution.observations(),
            Self::QuasiMonteCarlo(execution) => execution.observations(),
            Self::CopulaMonteCarlo(execution) => execution.monte_carlo().observations(),
            Self::CopulaQuasiMonteCarlo(execution) => execution.observations(),
        }
    }

    pub(super) fn evaluations_attempted(&self) -> usize {
        match self {
            Self::SobolSensitivity(execution) => execution.report().evaluations_attempted,
            Self::MonteCarlo(execution) => execution.evaluations_attempted(),
            Self::QuasiMonteCarlo(execution) => execution.evaluations_attempted(),
            Self::CopulaMonteCarlo(execution) => execution.monte_carlo().evaluations_attempted(),
            Self::CopulaQuasiMonteCarlo(execution) => execution.evaluations_attempted(),
        }
    }

    pub(super) fn status(&self) -> UqStatus {
        match self {
            Self::SobolSensitivity(execution) => execution.report().status,
            Self::MonteCarlo(execution) => execution.report().status,
            Self::QuasiMonteCarlo(execution) => execution.report().status,
            Self::CopulaMonteCarlo(execution) => execution.report().status,
            Self::CopulaQuasiMonteCarlo(execution) => execution.report().status,
        }
    }

    pub(super) fn rejection_reason(&self) -> Option<String> {
        match self {
            Self::SobolSensitivity(execution) => execution.report().rejection_reason,
            Self::MonteCarlo(execution) => execution.report().rejection_reason,
            Self::QuasiMonteCarlo(execution) => execution.report().rejection_reason,
            Self::CopulaMonteCarlo(execution) => execution.report().rejection_reason,
            Self::CopulaQuasiMonteCarlo(execution) => execution.report().rejection_reason,
        }
    }

    pub(super) fn monte_carlo(&self) -> Option<&UqExecution> {
        match self {
            Self::MonteCarlo(execution) => Some(execution),
            // Its observations/threshold are physical QoIs. Only its input
            // plan is latent; never use the inner sampler or checkpoint here.
            Self::CopulaMonteCarlo(execution) => Some(execution.monte_carlo()),
            Self::QuasiMonteCarlo(_) | Self::CopulaQuasiMonteCarlo(_) | Self::SobolSensitivity(_) => None,
        }
    }

    pub(super) fn qmc_report(&self) -> Option<QmcReport> {
        match self {
            Self::MonteCarlo(_) | Self::CopulaMonteCarlo(_) | Self::SobolSensitivity(_) => None,
            Self::QuasiMonteCarlo(execution) => Some(execution.report()),
            Self::CopulaQuasiMonteCarlo(execution) => Some(execution.report()),
        }
    }

    pub(super) fn sensitivity_report(&self) -> Option<SobolReport> {
        match self {
            Self::SobolSensitivity(execution) => Some(execution.report()),
            _ => None,
        }
    }

    pub(super) fn compliance(&self, model: &Model) -> Result<Option<compliance::Assessment>> {
        match self {
            Self::MonteCarlo(execution) => compliance::assess(model, execution),
            // Dependence is within each vector; Monte Carlo vectors remain
            // iid, so raw physical pass indicators retain Bernoulli semantics.
            Self::CopulaMonteCarlo(execution) => compliance::assess(model, execution.monte_carlo()),
            Self::QuasiMonteCarlo(_) | Self::CopulaQuasiMonteCarlo(_) | Self::SobolSensitivity(_) => Ok(None),
        }
    }

    pub(super) fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F)
    where F: FnMut(&[f64]) -> std::result::Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |values| evaluate(values).map(Some));
    }

    pub(super) fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, cancelled: C, evaluate: F,
    )
    where F: FnMut(&[f64]) -> std::result::Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        match self {
            Self::SobolSensitivity(execution) => {
                execution.advance_interruptible(allowance, cancelled, evaluate);
            }
            Self::MonteCarlo(execution) => {
                execution.advance_interruptible(allowance, cancelled, evaluate);
            }
            Self::QuasiMonteCarlo(execution) => {
                execution.advance_interruptible(allowance, cancelled, evaluate);
            }
            Self::CopulaMonteCarlo(execution) => {
                execution.advance_interruptible(allowance, cancelled, evaluate);
            }
            Self::CopulaQuasiMonteCarlo(execution) => {
                execution.advance_interruptible(allowance, cancelled, evaluate);
            }
        }
    }
}
