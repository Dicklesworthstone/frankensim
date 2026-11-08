//! Explicit resource limits for the existing sparse direct owner.
use super::LdltError;

/// Resource admission for cancellable AMD/symbolic/numeric preparation.
/// Counts bound owner work and payload, not process RSS or elapsed time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DirectLimits {
    /// Maximum input entries plus rows, and maximum symbolic L entries.
    pub max_symbolic_entries: usize,
    /// Cumulative graph visits, sort comparisons (conservative estimates),
    /// and index work during one symbolic analysis.
    pub max_analysis_work: usize,
    /// Conservative scalar-work allowance during one numeric factorization.
    pub max_numeric_work: usize,
    /// Largest admitted dense frontal dimension, including delayed pivots.
    pub max_front_dim: usize,
}

impl DirectLimits {
    pub(super) const UNLIMITED: Self = Self {
        max_symbolic_entries: usize::MAX,
        max_analysis_work: usize::MAX,
        max_numeric_work: usize::MAX,
        max_front_dim: usize::MAX,
    };
}

/// Internal accounting shared by the unchanged numeric/graph algorithms.
pub(super) struct DirectWork<'a> {
    limits: DirectLimits,
    cancelled: &'a dyn Fn() -> bool,
    numeric: bool,
    spent: usize,
}

impl<'a> DirectWork<'a> {
    pub(super) fn new(
        limits: DirectLimits,
        cancelled: &'a dyn Fn() -> bool,
        numeric: bool,
    ) -> Self {
        Self {
            limits,
            cancelled,
            numeric,
            spent: 0,
        }
    }

    pub(super) fn poll(&self) -> Result<(), LdltError> {
        if (self.cancelled)() {
            Err(LdltError::Cancelled)
        } else {
            Ok(())
        }
    }

    pub(super) fn charge(&mut self, amount: usize) -> Result<(), LdltError> {
        let cap = if self.numeric {
            self.limits.max_numeric_work
        } else {
            self.limits.max_analysis_work
        };
        let next = self.spent.checked_add(amount).ok_or(LdltError::WorkLimit {
            stage: if self.numeric { "numeric" } else { "symbolic" },
            required: usize::MAX,
            cap,
        })?;
        if next > cap {
            return Err(LdltError::WorkLimit {
                stage: if self.numeric { "numeric" } else { "symbolic" },
                required: next,
                cap,
            });
        }
        // Poll at least every 4096 accounted scalar/index operations. A
        // separately admitted native stage still owns its internal latency.
        if (next >> 12) != (self.spent >> 12) {
            self.poll()?;
        }
        self.spent = next;
        Ok(())
    }

    pub(super) fn symbolic(&self, entries: usize) -> Result<(), LdltError> {
        if entries > self.limits.max_symbolic_entries {
            return Err(LdltError::SymbolicBudget {
                entries,
                cap: self.limits.max_symbolic_entries,
            });
        }
        Ok(())
    }

    pub(super) fn front(&self, dimension: usize) -> Result<(), LdltError> {
        if dimension > self.limits.max_front_dim {
            return Err(LdltError::FrontBudget {
                dimension,
                cap: self.limits.max_front_dim,
            });
        }
        Ok(())
    }

    pub(super) fn bounded(&self) -> bool {
        self.limits != DirectLimits::UNLIMITED
    }

    pub(super) fn sort_work(length: usize) -> usize {
        length.saturating_mul(usize::BITS as usize - length.leading_zeros() as usize + 1)
    }
}
