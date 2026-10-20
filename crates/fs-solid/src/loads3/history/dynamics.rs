//! Feed existing typed load histories into an actual inertial trajectory.
//! Sampling remains in `TetLoadHistory3`; the sparse Newmark run owns time and
//! displacement/velocity/acceleration. These adapters never call a static solve.
use super::*;
use super::supports::SupportMotion3;
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

    /// Initialize driven supports from the SAME physical-time sample as every
    /// distributed force. Initial displacement and velocity are full node-major
    /// xyz vectors in world coordinates. Each restrained displacement must
    /// equal its sampled history value; omitted restrained components are zero.
    /// A disagreement refuses instead of silently overwriting an initial state.
    ///
    /// Initial restrained velocity and acceleration are explicit SI inputs,
    /// never derivatives guessed from two displacement samples. The full
    /// `prescribed_acceleration_m_s2` array must be zero on free components,
    /// whose acceleration is solved from the coupled inertial equations.
    /// Subsequent support velocities/accelerations use the prescribed-support
    /// Newmark recurrence through `advance_history_with_supports`.
    /// Thermal state, materials and reference constraint axes remain fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient_with_support_history<'g, 'm, 'p, 's>(
        &self,
        geometry: &'g TetLoadGeometry3<'m>,
        problem: &'p TetLinearElasticProblem<'s>,
        thermal: Option<TetThermalStrainField<'_>>,
        supports: &[SupportMotion3<'_>],
        time_s: f64,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        prescribed_acceleration_m_s2: &[[f64; 3]],
        config: TransientConfig3,
        cx: &Cx<'_>,
    ) -> Result<TetTransientRun3<'g, 'm, 'p, 's>, TransientHistoryError3> {
        cx.checkpoint().map_err(|_| LoadHistoryError3::Cancelled)?;
        if displacement_m.len() != problem.nodes_m.len() {
            return Err(TransientError3::Invalid(
                "initial support history needs one displacement vector per node",
            ).into());
        }
        let initial = self.at_with_supports(
            geometry, problem, supports, time_s, config.maximum_displacement_m, cx,
        )?;
        let prescribed = initial.prescribed_dofs_m();
        for &dof in problem.fixed_dofs {
            cx.checkpoint().map_err(|_| LoadHistoryError3::Cancelled)?;
            let expected = prescribed.binary_search_by_key(&dof, |&(d, _)| d)
                .map(|i| prescribed[i].1).unwrap_or(0.0);
            if displacement_m[dof / 3][dof % 3] != expected {
                return Err(TransientError3::Invalid(
                    "initial restrained displacement differs from its support history",
                ).into());
            }
        }
        Ok(initial.loads.case.begin_transient_with_supports(
            problem, thermal, time_s, displacement_m, velocity_m_s,
            prescribed_acceleration_m_s2, config, cx,
        )?)
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
    /// explicitly resolved forcing intervals. Driven supports use the explicit
    /// `advance_history_with_supports` path instead of this force-only adapter.
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

    /// Sample typed displacement and load histories together at the accepted
    /// run clock plus `h_s`, then execute ONE coupled inertial step. Declared
    /// frames use physical time even when a value table clamps or wraps.
    /// Restrained axes stay world Cartesian; rotating a payload rotates the
    /// prescribed displacement vector, not the constraint subspace.
    ///
    /// Missing fixed-component rows prescribe zero displacement. Support
    /// velocity/acceleration follow the discrete Newmark recurrence, so callers
    /// must resolve kinks and jumps with an appropriate step schedule. A value
    /// table is not a hidden analytic derivative or impulse model. There is no
    /// automatic event detection or temporal-error certificate.
    ///
    /// Failed sampling, invalid support assignments, cancellation and failed
    /// dynamic acceptance leave the previous load, support motion and clock
    /// intact. Reuse the same endpoint after correcting the failure.
    pub fn advance_history_with_supports(
        &mut self,
        history: &TetLoadHistory3<'_>,
        supports: &[SupportMotion3<'_>],
        h_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TransientStep3, TransientHistoryError3> {
        cx.checkpoint().map_err(|_| LoadHistoryError3::Cancelled)?;
        let time = self.state().time_s + h_s;
        if !h_s.is_finite() || h_s <= 0.0 || !time.is_finite() || time <= self.state().time_s {
            return Err(TransientError3::Invalid(
                "support history step must reach a finite later endpoint",
            ).into());
        }
        let endpoint = history.at_with_supports(
            self.geometry, self.problem, supports, time, self.maximum_displacement_m(), cx,
        )?;
        Ok(self.advance_with_supports(
            &endpoint.loads.case, endpoint.prescribed_dofs_m(), h_s, cx,
        )?)
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

impl crate::loads3::transient::adjoint::support::trajectory::RecordedPrescribedTransient3<'_, '_, '_, '_> {
    /// Sample typed forces and prescribed displacements together, then record
    /// one accepted coupled step. Reverse results address integrated forces and
    /// world-component support values at these exact physical times. Authored
    /// interpolation, frame and initial-velocity control derivatives are chained
    /// by the caller; this adapter does not guess them or discover events.
    pub fn advance_history_with_supports(
        &mut self, history: &TetLoadHistory3<'_>, supports: &[SupportMotion3<'_>],
        h_s: f64, cx: &Cx<'_>,
    ) -> Result<TransientStep3, crate::loads3::transient::adjoint::TransientAdjointError3> {
        self.admit_next_step(cx)?;
        let start = self.run().state().time_s;
        let time = start + h_s;
        if !h_s.is_finite() || h_s <= 0.0 || !time.is_finite() || time <= start {
            return Err(TransientError3::Invalid("recorded support history needs a finite later endpoint").into());
        }
        let endpoint = history.at_with_supports(
            self.run().geometry, self.run().problem, supports, time,
            self.run().maximum_displacement_m(), cx,
        )?;
        self.advance_with_supports(&endpoint.loads.case, endpoint.prescribed_dofs_m(), h_s, cx)
    }
}
