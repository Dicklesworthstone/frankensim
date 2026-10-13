//! External measurement acquisition, separate from numerical publication.
//!
//! Providers own real I/O, its timeouts, calibration and physical costs. They
//! receive only a selected request and an explicit execution context, never
//! hidden truth. A successful read is retained even if the post-I/O deadline or
//! cancellation check fails. Retrying numerical application does not read again.

use fs_exec::{AdmittedBudget, Cx};
use crate::OedError;
use super::measurement::{MeasurementError, MeasurementOutcome, MeasurementRequest, MeasurementStudy, MeasurementUpdate};

/// A real acquisition source. Implementations must honor the supplied context
/// and bound blocking I/O themselves: a synchronous callback cannot be forcibly
/// cancelled here. Use request.identity() as an external idempotency key where
/// supported. The framework never automatically retries a physical acquisition.
pub trait MeasurementProvider {
    type Error;
    fn read(&mut self, request: &MeasurementRequest, cx: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error>;
}

/// Owned acquired data that can be retried against its unchanged source study.
#[derive(Debug, Clone, PartialEq)]
pub struct AcquiredMeasurement {
    request: MeasurementRequest,
    outcome: MeasurementOutcome,
}
impl AcquiredMeasurement {
    #[must_use]
    pub const fn request(&self) -> &MeasurementRequest { &self.request }
    #[must_use]
    pub const fn outcome(&self) -> &MeasurementOutcome { &self.outcome }
    /// No provider is called. Failed application leaves both data and study
    /// available for retry; successful application makes the old request stale.
    pub fn apply<'a>(&self, study: &'a mut MeasurementStudy, cx: &Cx<'_>) -> Result<&'a MeasurementUpdate, MeasurementError> {
        study.assimilate(&self.request, &self.outcome, cx)
    }
    #[must_use]
    pub fn into_parts(self) -> (MeasurementRequest, MeasurementOutcome) { (self.request, self.outcome) }

    /// Encode a completed acquisition BEFORE applying it. Store these bytes
    /// durably with the source study checkpoint; no provider or I/O is called.
    pub fn to_checkpoint(&self, cx: &Cx<'_>) -> Result<Vec<u8>, MeasurementError> {
        super::measurement::checkpoint::acquisition_to_checkpoint(self, cx)
    }

    /// Recover an unapplied acquisition against its restored source study.
    /// Rechecks the current request; a stale or already-applied reading refuses.
    /// This never calls the provider, even when a previous application failed.
    pub fn from_checkpoint(study: &MeasurementStudy, bytes: &[u8], cx: &Cx<'_>) -> Result<Self, MeasurementError> {
        let (request, outcome) = super::measurement::checkpoint::acquisition_from_checkpoint(study, bytes, cx)?;
        Ok(Self { request, outcome })
    }
}

/// Distinguish failure before I/O, a provider error, and a real reading retained
/// after a numerical/protocol refusal. None of these variants refunds I/O costs.
#[derive(Debug)]
pub enum AcquisitionError<E> {
    /// No provider call was made.
    Before(MeasurementError),
    /// The provider returned no admitted reading. Its own side effects may
    /// already have happened; callers must not assume a retry is free or safe.
    Provider { request: MeasurementRequest, source: E },
    /// The provider returned a reading. Keep it for inspection or retry instead
    /// of losing it or acquiring the same sample a second time.
    After { acquired: AcquiredMeasurement, source: MeasurementError },
}
impl<E: core::fmt::Display> core::fmt::Display for AcquisitionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Before(e) => write!(f, "acquisition not started: {e}"),
            Self::Provider { source, .. } => write!(f, "measurement provider failed: {source}"),
            Self::After { source, .. } => write!(f, "acquired reading retained after refusal: {source}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for AcquisitionError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Before(e) | Self::After { source: e, .. } => Some(e),
            Self::Provider { source, .. } => Some(source),
        }
    }
}

impl MeasurementStudy {
    /// Acquire a selected current request WITHOUT changing the numerical study.
    /// The returned data is applied separately through AcquiredMeasurement::apply.
    /// All pre-I/O checks precede the provider call; all post-I/O failures retain
    /// the returned data. Planning, acquisition overhead, provider work and
    /// numerical application have separate explicit compute scopes.
    pub fn acquire<P: MeasurementProvider>(&self, request: &MeasurementRequest,
        provider: &mut P, cx: &Cx<'_>) -> Result<AcquiredMeasurement, AcquisitionError<P::Error>> {
        self.check_current(request).map_err(AcquisitionError::Before)?;
        let mut budget = AdmittedBudget::admit_ambient(cx, 2)
            .map_err(|e| AcquisitionError::Before(OedError::BudgetRefused(e).into()))?;
        budget.checkpoint("measurement before I/O", cx)
            .map_err(|e| AcquisitionError::Before(OedError::BudgetRefused(e).into()))?;
        budget.charge_cost("measurement before I/O", 1)
            .map_err(|e| AcquisitionError::Before(OedError::BudgetRefused(e).into()))?;
        let retained_request = request.clone();
        let outcome = match provider.read(request, cx) {
            Ok(outcome) => outcome,
            Err(source) => return Err(AcquisitionError::Provider { request: retained_request, source }),
        };
        let acquired = AcquiredMeasurement { request: retained_request, outcome };
        if acquired.outcome.request_identity() != acquired.request.identity() {
            return Err(AcquisitionError::After { acquired, source: MeasurementError::RequestMismatch });
        }
        let checked = budget.checkpoint("measurement after I/O", cx)
            .and_then(|()| budget.charge_cost("measurement after I/O", 1));
        if let Err(error) = checked {
            return Err(AcquisitionError::After { acquired, source: OedError::BudgetRefused(error).into() });
        }
        Ok(acquired)
    }
}

/// Explicit seeded-noise provider for synthetic acquisition experiments.
pub mod simulation;

/// Persistent acquisition allowances and single-use, pre-I/O reservations.
pub mod budget;
