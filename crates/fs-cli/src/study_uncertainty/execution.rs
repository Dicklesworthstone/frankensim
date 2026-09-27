//! Dispatch only: the statistical owners retain sampler, checkpoint and
//! failure semantics. Both paths evaluate the same native child-run callback.

use std::fmt::Display;

use fs_blake3::ContentHash;
use fs_uq::{QmcConfig, QmcExecution, QmcReport, UqExecution, UqStatus};

use super::{Model, Result, compliance, fail, plan};

pub(super) enum Execution {
    MonteCarlo(UqExecution),
    QuasiMonteCarlo(QmcExecution),
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
        if let Some(layout) = layout(model) {
            QmcExecution::new(&plan, layout).map(Self::QuasiMonteCarlo)
                .map_err(|error| fail("cli-uncertainty-plan", error))
        } else {
            UqExecution::new(&plan).map(Self::MonteCarlo)
                .map_err(|error| fail("cli-uncertainty-plan", error))
        }
    }

    pub(super) fn restore(model: &Model, bytes: &[u8]) -> Result<Self> {
        let plan = plan(model);
        if let Some(layout) = layout(model) {
            QmcExecution::restore(&plan, layout, model.identity(), bytes)
                .map(Self::QuasiMonteCarlo)
                .map_err(|error| fail("cli-uncertainty-resume", error.to_string()))
        } else {
            UqExecution::restore(&plan, model.identity(), bytes).map(Self::MonteCarlo)
                .map_err(|error| fail("cli-uncertainty-resume", error.to_string()))
        }
    }

    pub(super) fn checkpoint(&self, identity: ContentHash) -> Result<Vec<u8>> {
        let bytes = match self {
            Self::MonteCarlo(execution) => execution.checkpoint(identity),
            Self::QuasiMonteCarlo(execution) => execution.checkpoint(identity),
        };
        bytes.map_err(|error| fail("cli-uncertainty-checkpoint", error.to_string()))
    }

    pub(super) fn observations(&self) -> &[f64] {
        match self {
            Self::MonteCarlo(execution) => execution.observations(),
            Self::QuasiMonteCarlo(execution) => execution.observations(),
        }
    }

    pub(super) fn evaluations_attempted(&self) -> usize {
        match self {
            Self::MonteCarlo(execution) => execution.evaluations_attempted(),
            Self::QuasiMonteCarlo(execution) => execution.evaluations_attempted(),
        }
    }

    pub(super) fn status(&self) -> UqStatus {
        match self {
            Self::MonteCarlo(execution) => execution.report().status,
            Self::QuasiMonteCarlo(execution) => execution.report().status,
        }
    }

    pub(super) fn rejection_reason(&self) -> Option<String> {
        match self {
            Self::MonteCarlo(execution) => execution.report().rejection_reason,
            Self::QuasiMonteCarlo(execution) => execution.report().rejection_reason,
        }
    }

    pub(super) fn monte_carlo(&self) -> Option<&UqExecution> {
        match self {
            Self::MonteCarlo(execution) => Some(execution),
            Self::QuasiMonteCarlo(_) => None,
        }
    }

    pub(super) fn qmc_report(&self) -> Option<QmcReport> {
        match self {
            Self::MonteCarlo(_) => None,
            Self::QuasiMonteCarlo(execution) => Some(execution.report()),
        }
    }

    pub(super) fn compliance(&self, model: &Model) -> Result<Option<compliance::Assessment>> {
        match self {
            Self::MonteCarlo(execution) => compliance::assess(model, execution),
            Self::QuasiMonteCarlo(_) => Ok(None),
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
            Self::MonteCarlo(execution) => {
                execution.advance_interruptible(allowance, cancelled, evaluate);
            }
            Self::QuasiMonteCarlo(execution) => {
                execution.advance_interruptible(allowance, cancelled, evaluate);
            }
        }
    }
}
