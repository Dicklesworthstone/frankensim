//! Center-aligned complementary PWM on the existing independent gate contract.
//!
//! Each leg's ideal high pulse has width `d * period` and is centered in the
//! period. A dead interval of width `dead_time_s` is placed symmetrically about
//! each ideal transition: the outgoing device turns OFF half a dead time early,
//! and the incoming device turns ON half a dead time late. Thus an interior
//! duty has actual high/low ON times `d*T-t_dead` and `(1-d)*T-t_dead`.
//! This convention preserves pulse centering; it is not turn-ON-only delay.
//!
//! Exactly zero/one duty holds the low/high device continuously. An interior
//! pulse erased by blanking, an unrepresentable edge, or nonfinite data refuses
//! rather than clipping the duty or silently suppressing a pulse. A prepared
//! period is cyclic. Different periods can be joined through `admit_after`
//! only when their boundary gate commands agree; transitions to/from a held
//! opposite gate require an explicitly authored commutation interval.
//!
//! This fixed-size, allocation-free scheduler supplies commands, not a control
//! or loss model. Device currents, dead-time distortion and heat are resolved
//! by the original commutation owner. No stability, switching-transition loss,
//! asynchronous command-update or gate-driver propagation-delay claim is made.

use super::{CommutationInterval, DeviceGateState, invalid, positive};
use crate::PhsError;
use fs_math::dd::Dd;

/// Two legs have at most eight distinct transitions and nine positive intervals.
pub const MAX_PWM_INTERVALS: usize = 9;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CenterAlignedPwm {
    pub period_s: f64,
    /// Ideal high-device fractions for the input and output legs, in [0, 1].
    pub high_fraction: [f64; 2],
    /// One full both-OFF interval at each transition, independently for each leg.
    pub dead_time_s: [f64; 2],
}

/// Actual integrated command durations; these are not conducting-current times.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PwmTiming {
    /// Input high, input low, output high, output low.
    pub channel_on_s: [f64; 4],
    pub both_off_s: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub struct PwmProgram {
    command: CenterAlignedPwm,
    intervals: [CommutationInterval; MAX_PWM_INTERVALS],
    count: usize,
    boundary: DeviceGateState,
    timing: PwmTiming,
}

#[derive(Clone, Copy)]
struct Edge {
    time_s: f64,
    order: usize,
    leg: usize,
    on: [bool; 2],
}

impl CenterAlignedPwm {
    /// Bounded preparation: eight edge candidates, one fixed-size sort, no heap.
    pub fn prepare(self) -> Result<PwmProgram, PhsError> {
        if !positive(self.period_s)
            || self
                .high_fraction
                .iter()
                .any(|d| !d.is_finite() || !(0.0..=1.0).contains(d))
            || self
                .dead_time_s
                .iter()
                .any(|dt| !dt.is_finite() || *dt < 0.0)
        {
            return Err(invalid(
                "PWM requires a positive finite period, unit-interval duties and nonnegative finite dead times",
            ));
        }
        let mut edges = [Edge {
            time_s: 0.0,
            order: 0,
            leg: 0,
            on: [false; 2],
        }; 8];
        let mut edge_count = 0;
        let mut gates = DeviceGateState { on: [false; 4] };
        let period = Dd::from_f64(self.period_s);
        for leg in 0..2 {
            let duty = self.high_fraction[leg];
            gates.on[2 * leg] = duty == 1.0;
            gates.on[2 * leg + 1] = duty != 1.0;
            if duty == 0.0 || duty == 1.0 {
                continue;
            }
            let high = Dd::from_f64(self.period_s) * Dd::from_f64(duty);
            let low = period - high;
            let dead = Dd::from_f64(self.dead_time_s[leg]);
            if !positive((high - dead).to_f64()) || !positive((low - dead).to_f64()) {
                return Err(invalid(
                    "PWM dead time erases an interior high or low pulse",
                ));
            }
            let half = Dd::from_f64(0.5);
            let times = [
                ((low - dead) * half).to_f64(),
                ((low + dead) * half).to_f64(),
                (period - (low + dead) * half).to_f64(),
                (period - (low - dead) * half).to_f64(),
            ];
            if !positive(times[0])
                || times[0] > times[1]
                || times[1] >= times[2]
                || times[2] > times[3]
                || times[3] >= self.period_s
                || times.iter().any(|t| !t.is_finite())
                || (self.dead_time_s[leg] > 0.0 && (times[0] == times[1] || times[2] == times[3]))
            {
                return Err(invalid(
                    "PWM pulse or positive dead time is not representable on the period clock",
                ));
            }
            let states = [[false, false], [true, false], [false, false], [false, true]];
            for edge in 0..4 {
                edges[edge_count] = Edge {
                    time_s: times[edge],
                    order: 4 * leg + edge,
                    leg,
                    on: states[edge],
                };
                edge_count += 1;
            }
        }
        edges[..edge_count]
            .sort_unstable_by(|a, b| a.time_s.total_cmp(&b.time_s).then(a.order.cmp(&b.order)));
        let mut program = PwmProgram {
            command: self,
            intervals: [CommutationInterval {
                gates,
                duration_s: 0.0,
            }; MAX_PWM_INTERVALS],
            count: 0,
            boundary: gates,
            timing: PwmTiming {
                channel_on_s: [0.0; 4],
                both_off_s: [0.0; 2],
            },
        };
        let mut elapsed = Dd::ZERO;
        let mut i = 0;
        while i < edge_count {
            let edge_time = edges[i].time_s;
            program.append_to(edge_time, gates, &mut elapsed)?;
            // Apply every coincident edge before publishing the next interval.
            // Stable role order handles zero-dead-time OFF/ON pairs atomically.
            while i < edge_count && edges[i].time_s == edge_time {
                let edge = edges[i];
                gates.on[2 * edge.leg] = edge.on[0];
                gates.on[2 * edge.leg + 1] = edge.on[1];
                i += 1;
            }
        }
        program.append_to(self.period_s, gates, &mut elapsed)?;
        if gates != program.boundary {
            return Err(invalid(
                "PWM program does not close its periodic gate state",
            ));
        }
        let mut on = [Dd::ZERO; 4];
        let mut blank = [Dd::ZERO; 2];
        for interval in program.intervals() {
            let dt = Dd::from_f64(interval.duration_s);
            for (device, sum) in on.iter_mut().enumerate() {
                if interval.gates.on[device] {
                    *sum = *sum + dt;
                }
            }
            for (leg, sum) in blank.iter_mut().enumerate() {
                if !interval.gates.on[2 * leg] && !interval.gates.on[2 * leg + 1] {
                    *sum = *sum + dt;
                }
            }
        }
        program.timing = PwmTiming {
            channel_on_s: on.map(Dd::to_f64),
            both_off_s: blank.map(Dd::to_f64),
        };
        Ok(program)
    }
}

impl PwmProgram {
    fn append_to(
        &mut self,
        end_s: f64,
        gates: DeviceGateState,
        elapsed: &mut Dd,
    ) -> Result<(), PhsError> {
        let duration_s = (Dd::from_f64(end_s) - *elapsed).to_f64();
        if !positive(duration_s) || self.count >= MAX_PWM_INTERVALS {
            return Err(invalid(
                "PWM interval is nonpositive or exceeds its fixed edge capacity",
            ));
        }
        let next = *elapsed + Dd::from_f64(duration_s);
        if next.to_f64() != end_s {
            return Err(invalid(
                "PWM interval cannot retain its represented gate edge",
            ));
        }
        self.intervals[self.count] = CommutationInterval { gates, duration_s };
        self.count += 1;
        *elapsed = next;
        Ok(())
    }

    pub fn intervals(&self) -> &[CommutationInterval] {
        &self.intervals[..self.count]
    }
    pub fn command(&self) -> CenterAlignedPwm {
        self.command
    }
    pub fn period_s(&self) -> f64 {
        self.command.period_s
    }
    pub fn timing(&self) -> PwmTiming {
        self.timing
    }
    /// The same commands hold at both ends of this prepared periodic program.
    pub fn boundary_gates(&self) -> DeviceGateState {
        self.boundary
    }
    /// Admit direct concatenation without an instantaneous boundary gate change.
    /// This checks command continuity only; the native solver owns state/support.
    pub fn admit_after(&self, previous: &Self) -> Result<(), PhsError> {
        if self.boundary != previous.boundary {
            return Err(invalid(
                "PWM boundary gate change requires an explicit dead-time transition",
            ));
        }
        Ok(())
    }
}
