//! Passive local midpoint elimination of the original polarization/current pair.
//!
//! With h=dt/2 and d=1+h*gamma+(h*omega0)^2, midpoint current is
//! `Jm=(J0+h*omega0^2*(delta*Em-P0))/d`. Thus `P1=P0+dt*Jm`,
//! `J1=2*Jm-J0`, and the coefficient of E1 in P1 is the nonnegative
//! `delta*(h*omega0)^2/d`. This is the original two-state midpoint equation,
//! not an exponential propagator or a fitted relaxation law.

use super::{LorentzError, LorentzPole, LorentzState, Result, finite};
use fs_math::qd::Qd;

/// One fixed pole/interval condensation. All coefficients are prepared once;
/// local histories and recovery allocate nothing and retain original SI states.
#[derive(Clone, Copy, Debug)]
pub struct PreparedLorentzMidpoint {
    pole: LorentzPole,
    duration: Scaled,
    inverse_denominator: Scaled,
    current_gap_gain: Scaled,
    permittivity: Scaled,
    inverse_increment: Scaled,
    inverse_inertia: Scaled,
    heat_gain: Scaled,
}

/// One pole's recovered midpoint endpoint and independent energy densities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LorentzMidpointPoleStep<const N: usize> {
    pub state: LorentzState<N>,
    pub stored_before_j_m3: f64,
    pub stored_after_j_m3: f64,
    /// Quadratic change on actual published endpoints, not rounded-total subtraction.
    pub stored_change_j_m3: f64,
    /// Original midpoint electric work `Em dot (P1-P0)`.
    pub supplied_work_j_m3: f64,
    /// Independent `dt*gamma*|Jm|^2/(delta*omega0^2)`, before endpoint rounding.
    pub heat_j_m3: f64,
    /// Storage change plus heat minus work; includes published-state rounding.
    pub energy_residual_j_m3: f64,
}

impl LorentzPole {
    /// Condense this passive pole for one positive finite midpoint interval.
    /// No timestep stability cutoff is imposed. Phase/decay accuracy still
    /// requires resolving the pole and field; midpoint is second order and
    /// is not L-stable. Use the same preparation for history and recovery.
    pub fn prepare_midpoint(self, duration_s: f64) -> Result<PreparedLorentzMidpoint> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(LorentzError::Invalid(
                "positive finite Lorentz midpoint interval required",
            ));
        }
        let duration = Scaled::from(duration_s);
        let half = duration.mul(Scaled::from(0.5));
        let omega = Scaled::from(self.omega);
        let damping = Scaled::from(self.damping);
        let increment = Scaled::from(self.increment);
        let a = half.mul(omega);
        let denominator = Scaled::from(1.0).add(half.mul(damping)).add(a.mul(a));
        let inverse_denominator = denominator.recip();
        let current_gap_gain = half.mul(omega).mul(omega).mul(inverse_denominator);
        let inverse_increment = increment.recip();
        let inverse_inertia = increment.mul(omega).mul(omega).recip();
        Ok(PreparedLorentzMidpoint {
            pole: self,
            duration,
            inverse_denominator,
            current_gap_gain,
            permittivity: half.mul(current_gap_gain).mul(increment),
            inverse_increment,
            inverse_inertia,
            heat_gain: duration.mul(damping).mul(inverse_inertia),
        })
    }
}

impl PreparedLorentzMidpoint {
    /// Nonnegative endpoint-field coefficient in `P1 = gain*E1 + history` [F/m].
    pub fn permittivity_increment_f_m(&self) -> Result<f64> {
        self.permittivity.physical()
    }

    /// Original history in the condensed endpoint polarization [C/m2].
    pub fn history_c_m2<const N: usize>(
        &self,
        electric_before_v_m: [f64; N],
        before: LorentzState<N>,
    ) -> Result<[f64; N]> {
        admit(before, electric_before_v_m, [0.0; N])?;
        let mut history = [0.0; N];
        for i in 0..N {
            let current = self.midpoint_current(before, i, electric_before_v_m[i], 0.0);
            history[i] = Scaled::from(before.polarization_c_m2[i])
                .add(self.duration.mul(current))
                .physical()?;
        }
        Ok(history)
    }

    /// Recover P/J after the coupled field solve without changing the input.
    /// Heat is independently evaluated from the unrounded midpoint current;
    /// it is never assigned a work residual. State publication contributes its
    /// own reported balance residual. Arbitrary large dt is passive, but does
    /// not imply resolved physical oscillation or accurate material heating.
    pub fn advance<const N: usize>(
        &self,
        before: LorentzState<N>,
        electric_before_v_m: [f64; N],
        electric_after_v_m: [f64; N],
    ) -> Result<LorentzMidpointPoleStep<N>> {
        admit(before, electric_before_v_m, electric_after_v_m)?;
        let stored_before = self.pole.stored_energy_j_m3(before)?;
        let mut state = LorentzState::default();
        let mut change = Scaled::ZERO;
        let mut work = Scaled::ZERO;
        let mut heat = Scaled::ZERO;
        for i in 0..N {
            let p0 = Scaled::from(before.polarization_c_m2[i]);
            let j0 = Scaled::from(before.current_a_m2[i]);
            let jm =
                self.midpoint_current(before, i, electric_before_v_m[i], electric_after_v_m[i]);
            state.polarization_c_m2[i] = p0.add(self.duration.mul(jm)).physical()?;
            state.current_a_m2[i] = jm.mul(Scaled::from(2.0)).add(j0.neg()).physical()?;
            let p1 = Scaled::from(state.polarization_c_m2[i]);
            let j1 = Scaled::from(state.current_a_m2[i]);
            let dp = p1.add(p0.neg());
            let dj = j1.add(j0.neg());
            change = change
                .add(
                    p0.add(p1)
                        .mul(Scaled::from(0.5))
                        .mul(dp)
                        .mul(self.inverse_increment),
                )
                .add(
                    j0.add(j1)
                        .mul(Scaled::from(0.5))
                        .mul(dj)
                        .mul(self.inverse_inertia),
                );
            work = work.add(mean(electric_before_v_m[i], electric_after_v_m[i]).mul(dp));
            heat = heat.add(jm.mul(jm).mul(self.heat_gain));
        }
        let heat_j_m3 = heat.physical()?;
        if heat_j_m3 < 0.0 {
            return Err(LorentzError::Invalid(
                "Lorentz midpoint heat lost nonnegativity",
            ));
        }
        Ok(LorentzMidpointPoleStep {
            state,
            stored_before_j_m3: stored_before,
            stored_after_j_m3: self.pole.stored_energy_j_m3(state)?,
            stored_change_j_m3: change.physical()?,
            supplied_work_j_m3: work.physical()?,
            heat_j_m3,
            energy_residual_j_m3: change.add(heat).add(work.neg()).physical()?,
        })
    }

    fn midpoint_current<const N: usize>(
        &self,
        before: LorentzState<N>,
        i: usize,
        field_before: f64,
        field_after: f64,
    ) -> Scaled {
        // Cancel the original restoring gap before adding the current memory.
        let gap = Scaled::from(self.pole.increment)
            .mul(mean(field_before, field_after))
            .add(Scaled::from(before.polarization_c_m2[i]).neg());
        self.inverse_denominator
            .mul(Scaled::from(before.current_a_m2[i]))
            .add(self.current_gap_gain.mul(gap))
    }
}

fn admit<const N: usize>(before: LorentzState<N>, e0: [f64; N], e1: [f64; N]) -> Result<()> {
    if !(1..=3).contains(&N)
        || before
            .polarization_c_m2
            .iter()
            .chain(&before.current_a_m2)
            .chain(&e0)
            .chain(&e1)
            .any(|v| !v.is_finite())
    {
        Err(LorentzError::Invalid(
            "finite 1D/2D/3D Lorentz state and field vectors required",
        ))
    } else {
        Ok(())
    }
}

// Source-bank instantaneous field storage/work use the same range-safe
// products; the source wrapper owns epsilon_infinity and all state admission.
pub(crate) fn field_energy<const N: usize>(epsilon: f64, field: [f64; N]) -> Result<f64> {
    let mut sum = Scaled::ZERO;
    for e in field {
        sum = sum.add(
            Scaled::from(epsilon)
                .mul(Scaled::from(e))
                .mul(Scaled::from(e))
                .mul(Scaled::from(0.5)),
        );
    }
    sum.physical()
}

pub(crate) fn field_increment<const N: usize>(
    epsilon: f64,
    before: [f64; N],
    after: [f64; N],
) -> Result<f64> {
    let mut sum = Scaled::ZERO;
    for i in 0..N {
        let delta = Scaled::from(after[i]).add(Scaled::from(before[i]).neg());
        sum = sum.add(
            Scaled::from(epsilon)
                .mul(mean(before[i], after[i]))
                .mul(delta),
        );
    }
    sum.physical()
}

fn mean(a: f64, b: f64) -> Scaled {
    Scaled::from(a).add(Scaled::from(b)).mul(Scaled::from(0.5))
}

// Local dimensional products keep an explicit binary exponent until their
// physical result is published. The existing quad-double owner retains a
// nearly balanced restoring gap and midpoint endpoint increments. This avoids
// squaring dimensional optical rates or forming an overflowing dt*omega.
// It extends range, not floating precision; no interval enclosure is claimed.
#[derive(Clone, Copy, Debug)]
struct Scaled {
    value: Qd,
    exponent: i32,
}

impl Scaled {
    const ZERO: Self = Self {
        value: Qd::ZERO,
        exponent: 0,
    };
    fn from(value: f64) -> Self {
        Self::normalized(Qd::from_f64(value), 0)
    }
    fn normalized(value: Qd, exponent: i32) -> Self {
        if value.is_zero() {
            return Self::ZERO;
        }
        let lead = value.to_f64().abs();
        let (normal, correction) = if lead < f64::MIN_POSITIVE {
            (lead * 18_014_398_509_481_984.0, -54)
        } else {
            (lead, 0)
        };
        let shift = ((normal.to_bits() >> 52) & 0x7ff) as i32 - 1023 + correction;
        Self {
            value: value.scale_power_of_two(-shift),
            exponent: exponent + shift,
        }
    }
    fn neg(self) -> Self {
        Self {
            value: -self.value,
            ..self
        }
    }
    fn add(self, other: Self) -> Self {
        if self.value.is_zero() {
            return other;
        }
        if other.value.is_zero() {
            return self;
        }
        let exponent = self.exponent.max(other.exponent);
        Self::normalized(
            self.value.scale_power_of_two(self.exponent - exponent)
                + other.value.scale_power_of_two(other.exponent - exponent),
            exponent,
        )
    }
    fn mul(self, other: Self) -> Self {
        if self.value.is_zero() || other.value.is_zero() {
            return Self::ZERO;
        }
        Self::normalized(self.value * other.value, self.exponent + other.exponent)
    }
    fn recip(self) -> Self {
        Self::normalized(Qd::ONE / self.value, -self.exponent)
    }
    fn physical(self) -> Result<f64> {
        finite(self.value.scale_power_of_two(self.exponent).to_f64())
    }
}

#[cfg(test)]
#[path = "midpoint_tests.rs"]
mod tests;
