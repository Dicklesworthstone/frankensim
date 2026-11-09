//! Shared bounded gate/refinement plan; no electrical state or solver lives here.
use super::*;
pub(super) struct GateStep {
    pub(super) gates: DeviceGateState,
    pub(super) duration_s: f64,
    pub(super) legs: [HalfBridgeConduction; 2],
}
pub(super) struct GatePlan {
    pub(super) intervals: Vec<CommutationInterval>,
    pub(super) indices: Vec<[usize; 2]>,
    pub(super) caches: Vec<GateStep>,
    pub(super) substeps: usize,
    pub(super) integration_steps: usize,
    pub(super) period_s: f64,
    pub(super) electrical: PowerLoadControls,
    pub(super) controls: SwitchCycleControls,
}
impl GatePlan {
    pub(super) fn new(
        devices: [SwitchDevice; 4],
        intervals: &[CommutationInterval],
        substeps: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(cancel)?;
        if !(1..=MAX_SWITCH_INTERVALS).contains(&controls.max_intervals)
            || intervals.is_empty()
            || intervals.len() > controls.max_intervals
            || !(1..=MAX_POWER_SUBSTEPS_PER_INTERVAL).contains(&substeps)
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
            || !electrical.absolute_energy_tolerance_j.is_finite()
            || electrical.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&electrical.relative_energy_tolerance)
        {
            return Err(invalid("capacitive commutation needs bounded intervals/refinement and finite energy controls").into());
        }
        electrical.solve.validate()?;
        let integration_steps = intervals
            .len()
            .checked_mul(substeps)
            .ok_or_else(|| invalid("capacitive integration count overflow"))?;
        integration_steps
            .checked_mul(electrical.solve.max_iterations)
            .and_then(|n| n.checked_mul(usize::from(electrical.backtracking_halvings) + 1))
            .ok_or_else(|| invalid("capacitive nonlinear work bound overflow"))?;
        let mut period = Dd::ZERO;
        let mut caches: Vec<GateStep> = Vec::new();
        let mut indices = Vec::with_capacity(intervals.len());
        for interval in intervals {
            poll(cancel)?;
            if !positive(interval.duration_s) {
                return Err(invalid("positive finite capacitive gate duration required").into());
            }
            period = period + Dd::from_f64(interval.duration_s);
            let regular = interval.duration_s / substeps as f64;
            let last = number(
                Dd::from_f64(interval.duration_s)
                    - Dd::from_f64(regular) * Dd::from_f64((substeps - 1) as f64),
            )?;
            if !positive(regular) || !positive(last) {
                return Err(invalid("capacitive subdivision duration is unrepresentable").into());
            }
            let mut represented = Dd::ZERO;
            for j in 0..substeps {
                poll(cancel)?;
                represented =
                    represented + Dd::from_f64(if j + 1 == substeps { last } else { regular });
            }
            if number(represented)?.to_bits() != interval.duration_s.to_bits() {
                return Err(invalid(
                    "capacitive subdivisions must preserve each authored gate edge",
                )
                .into());
            }
            let mut pair = [0; 2];
            for (slot, duration_s) in [regular, last].into_iter().enumerate() {
                poll(cancel)?;
                pair[slot] = if let Some(i) = caches.iter().position(|c| {
                    c.gates == interval.gates && c.duration_s.to_bits() == duration_s.to_bits()
                }) {
                    i
                } else {
                    let on = interval.gates.on;
                    let legs = [
                        HalfBridgeConduction::new(devices[0], devices[1], [on[0], on[1]])?,
                        HalfBridgeConduction::new(devices[2], devices[3], [on[2], on[3]])?,
                    ];
                    let i = caches.len();
                    caches.push(GateStep {
                        gates: interval.gates,
                        duration_s,
                        legs,
                    });
                    i
                };
            }
            indices.push(pair);
        }
        let period_s = number(period)?;
        if !positive(period_s) {
            return Err(invalid("capacitive cycle duration is unrepresentable").into());
        }
        poll(cancel)?;
        Ok(Self {
            intervals: intervals.to_vec(),
            indices,
            caches,
            substeps,
            integration_steps,
            period_s,
            electrical,
            controls,
        })
    }
    pub(super) fn step(&self, k: usize, j: usize) -> &GateStep {
        &self.caches[self.indices[k][usize::from(j + 1 == self.substeps)]]
    }
    pub(super) fn workspace(&self, system: &PortHamiltonian) -> Result<StepWorkspace, PhsError> {
        let mut workspace = StepWorkspace::new(system)?;
        workspace.set_solve_budget(self.electrical.solve)?;
        workspace.set_backtracking_limit(self.electrical.backtracking_halvings)?;
        Ok(workspace)
    }
}
