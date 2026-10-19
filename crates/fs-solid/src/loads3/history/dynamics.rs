//! Feed existing typed load histories into an actual inertial trajectory.
//! Sampling remains in `TetLoadHistory3`; the sparse Newmark run owns time and
//! displacement/velocity/acceleration. These adapters never call a static solve.
use super::*;
use crate::loads3::transient::{TetTransientRun3, TransientConfig3, TransientError3, TransientStep3};

#[derive(Debug, Clone, PartialEq)]
pub enum TransientHistoryError3 {
    History(LoadHistoryError3),
    Dynamics(TransientError3),
}
impl std::fmt::Display for TransientHistoryError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "scheduled tetrahedral dynamics refused: {self:?}")
    }
}
impl std::error::Error for TransientHistoryError3 {}
impl From<LoadHistoryError3> for TransientHistoryError3 {
    fn from(e: LoadHistoryError3) -> Self { Self::History(e) }
}
impl From<TransientError3> for TransientHistoryError3 {
    fn from(e: TransientError3) -> Self { Self::Dynamics(e) }
}
impl TetLoadHistory3<'_> {
    /// Sample the initial source at its physical time and use that same source
    /// to initialize M a0=f0-C v0-K u0. Density for acceleration histories is the
    /// exact borrowed solve material, not an independent mass approximation.
    /// Thermal strain, materials, mesh and zero-valued supports remain frozen.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient<'g, 'm, 'p, 's>(&self, geometry: &'g TetLoadGeometry3<'m>,
        problem: &'p TetLinearElasticProblem<'s>, thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64, displacement_m: &[[f64; 3]], velocity_m_s: &[[f64; 3]],
        config: TransientConfig3, cx: &Cx<'_>)
        -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientHistoryError3> {
        let initial = self.at(geometry, problem, time_s, cx)?;
        Ok(initial.case.begin_transient(problem, thermal, time_s, displacement_m,
            velocity_m_s, config, cx)?)
    }
}
impl TetTransientRun3<'_, '_, '_, '_> {
    /// Sample every force at the Newmark endpoint t+h, preserving the physical
    /// frame clock even when a value table clamps/wraps. Reuse prepared K/M;
    /// only load integration and one inertial step are performed. A failed source
    /// or step leaves the accepted motion, clock and previous forcing unchanged.
    ///
    /// Endpoint forces are integrated by Newmark's trapezoidal work rule.
    /// This does NOT detect events, resolve jumps or certify temporal error.
    /// Smooth ramps need a time-step ladder; discontinuous/short pulses need
    /// explicitly resolved forcing intervals. Existing support-motion histories
    /// are quasi-static and are not silently treated as inertial base excitation.
    pub fn advance_history(&mut self, history: &TetLoadHistory3<'_>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientStep3, TransientHistoryError3> {
        if !h_s.is_finite() || h_s <= 0.0 {
            return Err(TransientError3::Invalid("history step must be finite and positive").into());
        }
        let time = self.state().time_s + h_s;
        if !time.is_finite() || time <= self.state().time_s {
            return Err(TransientError3::Invalid("history endpoint must be finite and later").into());
        }
        let geometry = self.geometry; let problem = self.problem;
        let endpoint = history.at(geometry, problem, time, cx)?;
        Ok(self.advance(&endpoint.case, h_s, cx)?)
    }
}

/// Finite schedule indexing and exact left/right transient load events.
pub mod events;

impl crate::loads3::transient::adjoint::trajectory::RecordedTransient3<'_, '_, '_, '_> {
    /// Use the original typed source sampler, retaining its actual endpoint
    /// load for control differentiation. The reverse result is with respect to
    /// integrated endpoint forces; source parameters chain through load
    /// directions/interpolation weights supplied by their owner. This method
    /// does not discover events or differentiate source frame/time parameters.
    pub fn advance_history(&mut self, history: &TetLoadHistory3<'_>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientStep3, crate::loads3::transient::adjoint::TransientAdjointError3> {
        let start = self.run().state().time_s;
        let time = start + h_s;
        if !h_s.is_finite() || h_s <= 0.0 || !time.is_finite() || time <= start {
            return Err(TransientError3::Invalid("recorded history step must reach a finite later time").into());
        }
        let geometry = self.run().geometry; let problem = self.run().problem;
        let endpoint = history.at(geometry, problem, time, cx)?;
        self.advance(&endpoint.case, h_s, cx)
    }
}
