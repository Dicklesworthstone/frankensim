//! Explicit finite load jumps, with continuous displacement and velocity.
//! An impulse is a different physical operation and is NOT represented here.
use super::*;

/// Diagnostics for changing a finite force at the current accepted time.
/// No time, step, displacement, velocity or energy is changed by this operation.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadJump3 {
    pub time_s: f64,
    pub iterations: usize,
    pub force_residual_n: f64,
    /// Euclidean norm over free xyz acceleration components; not a time error.
    pub acceleration_change_l2_m_s2: f64,
}

/// One original Newmark interval to the event's LEFT load, followed by the
/// zero-duration switch to its RIGHT load. Both operations publish atomically.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientEventStep3 {
    pub interval: TransientStep3,
    pub jump: LoadJump3,
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Apply a finite load discontinuity at the CURRENT time. Solve
    /// M*a_plus = f_plus + f_thermal - C*v - K*u using the original mass solve.
    /// A finite jump changes acceleration, not momentum: q/v and both energy
    /// fields are preserved bit-for-bit. No impulse or artificial work is added.
    /// The full new load, including fixed-node forces, becomes the accepted
    /// response input and the next Newmark interval's starting load.
    ///
    /// No step slot is consumed because time does not advance. Work remains
    /// bounded by the existing linear-iteration allowance. Source/solve failure
    /// or cancellation publishes neither the new force nor acceleration.
    pub fn apply_load_jump(&mut self, load_after: &TetLoadCase3<'_, '_>, cx: &Cx<'_>)
        -> Result<LoadJump3, TransientError3> {
        check(cx)?;
        let force = self.reduce_force(load_after, cx)?;
        let ku = self.apply(0.0, 1.0, &self.state.displacement, cx)?;
        let cv = self.apply(self.config.mass_damping_per_s, self.config.stiffness_damping_s,
            &self.state.velocity, cx)?;
        let mut rhs = Vec::with_capacity(force.len());
        for i in 0..force.len() {
            if i % 256 == 0 { check(cx)?; }
            rhs.push(number(force[i] + self.thermal[i] - ku[i] - cv[i])?);
        }
        let (a, iterations, force_residual_n) = solve(&self.operator(1.0, 0.0, cx), &rhs, self.config)?;
        let mut change = 0.0_f64;
        for (i, (&next, &old)) in a.iter().zip(&self.state.acceleration).enumerate() {
            if i % 256 == 0 { check(cx)?; }
            change = number(change.hypot(number(next-old)?))?;
        }
        let mut full_force = Vec::with_capacity(load_after.nodal_forces_n().len());
        for &(d, f) in load_after.nodal_forces_n() { check(cx)?; full_force.push((d, f)); }
        let report = LoadJump3 { time_s: self.state.time_s, iterations, force_residual_n,
            acceleration_change_l2_m_s2: change };
        check(cx)?;
        self.state.acceleration = a;
        self.force = force;
        self.full_force = full_force;
        Ok(report)
    }

    /// Advance to an explicitly resolved event, then switch finite loads there.
    /// `load_before` and `load_after` are the physical LEFT and RIGHT limits at
    /// the SAME endpoint. Use the left load for the preceding interval; using
    /// the right load there smears a step backwards through the entire interval.
    /// Unchanged Newmark quadrature applies between the retained starting load
    /// and `load_before`. Any intervening discontinuities must be split by the
    /// caller. This does not discover knots, approximate left limits by epsilon,
    /// wrap frame clocks, or claim to detect unknown sub-step pulses.
    ///
    /// If either the interval or right-limit mass solve fails, neither is
    /// published. Forking shares the assembled operators and frozen thermal data.
    /// A successful event consumes exactly one ordinary advancement step slot.
    pub fn advance_to_load_jump(&mut self, load_before: &TetLoadCase3<'_, '_>,
        load_after: &TetLoadCase3<'_, '_>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientEventStep3, TransientError3> {
        check(cx)?;
        let mut candidate = self.clone();
        let interval = candidate.advance(load_before, h_s, cx)?;
        let jump = candidate.apply_load_jump(load_after, cx)?;
        check(cx)?;
        *self = candidate;
        Ok(TransientEventStep3 { interval, jump })
    }
}
