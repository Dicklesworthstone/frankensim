//! Event-aligned native thermal time grids shared by admission and execution.

use std::collections::BTreeSet;

use fs_scenario::Violation;

use crate::spec::{ConductionTransient, dims};

fn invalid(code: &'static str, what: impl Into<String>, fix: &str) -> Violation {
    Violation {
        code,
        what: what.into(),
        fix: fix.to_string(),
    }
}

fn work_limit() -> Violation {
    invalid(
        "project-conduction-transient-steps",
        "the event-aligned coarse grid and its nested half grid exceed the declared combined step cap",
        "increase max-steps to cover three times the coarse intervals, at most 10000; each workload switch consumes a time-grid boundary",
    )
}

impl ConductionTransient {
    /// Build coarse endpoints including every declared workload switch.
    ///
    /// Each segment between switches is divided into
    /// `ceil(segment_duration / max_step)` intervals. The nested fine grid
    /// bisects each returned interval, so their combined work is three times
    /// this vector's length. With no schedules, this retains the original
    /// uniform grid's endpoint arithmetic. At most `max_steps / 3` endpoints
    /// and distinct switch times are allocated, regardless of input size.
    ///
    /// # Errors
    /// Refuses invalid time units, nonpositive or nonfinite controls, invalid
    /// power intervals, an insufficient combined work cap, or a coarse/fine
    /// interval whose positive duration cannot be represented. Region ownership
    /// and source coverage are additionally checked by `ProjectSpec::validate`.
    pub fn coarse_step_ends_s(&self) -> Result<Vec<f64>, Violation> {
        let horizon = self.horizon.value;
        let max_step = self.max_step.value;
        if self.horizon.dims != dims::TIME
            || self.max_step.dims != dims::TIME
            || !(horizon.is_finite() && horizon > 0.0 && max_step.is_finite() && max_step > 0.0)
        {
            return Err(invalid(
                "project-conduction-transient-quantity",
                "horizon and max-step must be positive finite quantities in seconds",
                "declare the physical horizon and maximum coarse step in coherent SI time units",
            ));
        }
        if !(3..=10_000).contains(&self.max_steps) {
            return Err(work_limit());
        }
        let cap = (self.max_steps / 3) as usize;
        // Positive finite f64 bit patterns have the same ordering as their
        // values. This set unifies simultaneous switches without tolerances
        // that could erase a short but explicitly declared pulse.
        let mut switches = BTreeSet::new();
        switches.insert(horizon.to_bits());
        for schedule in &self.power_schedules {
            if schedule.steps.is_empty() {
                return Err(invalid(
                    "project-conduction-transient-power-step",
                    format!("power schedule for `{}` has no intervals", schedule.region),
                    "declare one or more strictly increasing interval endpoints ending at the horizon",
                ));
            }
            // One schedule's strict endpoints are distinct, hence each needs
            // at least one coarse interval before any union is constructed.
            if schedule.steps.len() > cap {
                return Err(work_limit());
            }
            let mut previous = 0.0;
            for step in &schedule.steps {
                let until = step.until.value;
                if step.until.dims != dims::TIME
                    || !until.is_finite()
                    || until <= previous
                    || until > horizon
                    || step.watts.dims != dims::POWER
                    || !step.watts.value.is_finite()
                    || step.watts.value < 0.0
                {
                    return Err(invalid(
                        "project-conduction-transient-power-step",
                        format!(
                            "power schedule for `{}` needs increasing endpoints within the horizon and finite nonnegative delivered watts",
                            schedule.region
                        ),
                        "use seconds for until and watts for power; every interval must have positive duration",
                    ));
                }
                if switches.len() == cap && !switches.contains(&until.to_bits()) {
                    return Err(work_limit());
                }
                switches.insert(until.to_bits());
                previous = until;
            }
            if previous != horizon {
                return Err(invalid(
                    "project-conduction-transient-power-horizon",
                    format!(
                        "power schedule for `{}` ends at {previous} s, not the declared {horizon} s horizon",
                        schedule.region
                    ),
                    "end every regional workload exactly at the horizon; no tail power is inferred",
                ));
            }
        }
        let mut endpoints = Vec::with_capacity(cap);
        let mut start = 0.0;
        for bits in switches {
            let end = f64::from_bits(bits);
            let duration = end - start;
            let count = (duration / max_step).ceil().max(1.0);
            if !count.is_finite() || count > (cap - endpoints.len()) as f64 {
                return Err(work_limit());
            }
            let count = count as usize;
            let mut previous = start;
            for ordinal in 1..=count {
                let time = if ordinal == count {
                    end
                } else {
                    start + duration * (ordinal as f64 / count as f64)
                };
                let midpoint = previous + 0.5 * (time - previous);
                if !(time > previous && time <= end && midpoint > previous && midpoint < time) {
                    return Err(invalid(
                        "project-conduction-transient-time-resolution",
                        "a requested coarse interval or its half step has no representably positive duration",
                        "choose representable switch times and a time step that supports a nested half grid; no pulse can be skipped",
                    ));
                }
                endpoints.push(time);
                previous = time;
            }
            start = end;
        }
        Ok(endpoints)
    }
}
