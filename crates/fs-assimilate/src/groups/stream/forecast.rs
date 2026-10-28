//! Advance the physical state while preserving learned calibration covariance.

use fs_blake3::DomainHasher;

use super::{SharedNoiseStream, number};
use crate::factor::forecast::{LinearDynamics, invalid};
use crate::{AssimError, Cx, WorkPlan, WorkProgress, canonicalize_zero};

impl SharedNoiseStream {
    /// Advance `x` through a complete discrete linear-Gaussian step while
    /// retaining the SAME calibration variables `b` and learned `Cov(x,b)`.
    /// The joint law is `[x_next,b] = diag(F,I)[x,b] + [offset,0] + [w,0]`.
    /// Process noise must be independent of BOTH the physical and nuisance
    /// priors; nuisance drift and cross-correlated process noise are not modeled.
    ///
    /// The first successful call anchors the existing state at `from_time`.
    /// Thereafter `from_time` must exactly equal `time()` and `to_time` must be
    /// finite and later. Times are caller-chosen physical units; the transition
    /// already represents the entire interval. It is never multiplied by dt.
    ///
    /// Prediction changes the history root, invalidating previously planned
    /// sensor acquisitions. Accepted observation IDs remain retained, so old
    /// readings cannot acquire a second update merely because time advanced.
    /// New readings passed to `assimilate` are assumed current; delayed data
    /// requires a separate smoothing model. Forecasts themselves are not
    /// replayed: a repeated interval fails the frontier check.
    ///
    /// All fallible work precedes publication of the posterior, history and
    /// time together. Budget/cancellation/arithmetic refusals preserve all three.
    /// Prediction and history preparation retain separate checked Cx scopes,
    /// matching the stream's existing update composition contract.
    pub fn forecast_physical(
        &mut self,
        dynamics: &LinearDynamics,
        from_time: f64,
        to_time: f64,
        cx: &Cx<'_>,
    ) -> Result<(), AssimError> {
        if dynamics.dimension() != self.physical_dim {
            return Err(AssimError::DimMismatch {
                state: self.physical_dim,
                operator: dynamics.dimension(),
            });
        }
        if !from_time.is_finite()
            || !to_time.is_finite()
            || to_time <= from_time
            || !(to_time - from_time).is_finite()
        {
            return Err(invalid(
                "forecast times must be finite, increasing, with finite duration",
            ));
        }
        if self.time.is_some_and(|current| current != from_time) {
            return Err(invalid(
                "forecast start differs from the retained time frontier",
            ));
        }
        let du = self.joint.dim() as u128;
        let mut progress = WorkProgress::new(
            cx,
            WorkPlan::checked(0, 0, 0, 0, 0, 0, 16 * du * du + 32 * du + 128)?,
        )?;
        progress.checkpoint("shared stream forecast admission")?;
        let joint = dynamics.predict_prefix(&self.joint, cx)?;
        let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.shared-stream-forecast.v1");
        hash.update(self.history.as_bytes());
        hash.update(dynamics.content_id().as_bytes());
        progress.hash_bytes("shared stream forecast identity", 64)?;
        number(&mut hash, from_time, &mut progress)?;
        number(&mut hash, to_time, &mut progress)?;
        // Bind the actual numerical result as well as its producing inputs.
        for (row, value) in joint.mean().iter().enumerate() {
            number(&mut hash, *value, &mut progress)?;
            number(
                &mut hash,
                joint.diag(row).expect("bounded factor row"),
                &mut progress,
            )?;
            for column in row + 1..joint.dim() {
                number(
                    &mut hash,
                    joint.upper_entry(row, column).expect("strict upper entry"),
                    &mut progress,
                )?;
            }
        }
        let history = hash.finalize();
        progress.checkpoint("shared stream forecast publication")?;
        self.joint = joint;
        self.history = history;
        self.time = Some(canonicalize_zero(to_time));
        Ok(())
    }
}
